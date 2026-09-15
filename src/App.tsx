/**
 * @license
 * SPDX-License-Identifier: Apache-2.0
 */

import { useState, FormEvent, useEffect, useRef } from "react";
import { motion, AnimatePresence } from "motion/react";
import { 
  Search, 
  AlertTriangle, 
  CheckCircle2, 
  ShieldAlert,
  RefreshCw,
  Info,
  Database,
  History,
  ArrowLeft,
  FileSpreadsheet
} from "lucide-react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle, CardFooter } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { predictInteraction, PredictionResult, InputType, CompoundInput, AnalysisError, PredictionMethod, FunctionalGroupReactivity } from "@/src/lib/gemini";
import { detectFunctionalGroupsDetailed } from "@/src/lib/reaction-engine";
import { ChemicalStructure } from "@/src/components/ChemicalStructure";
import { InteractionHeatmap } from "@/src/components/InteractionHeatmap";
import { CsvLogbook } from "@/src/components/CsvLogbook";
import { initRDKit, getMolecularDescriptors, computeStrainEnergy } from "@/src/lib/rdkit";
import { sanitizeData } from "@/src/lib/firestore-utils";
import { Plus, Trash2, AlertCircle, WifiOff, Clock, Lock } from "lucide-react";
import * as XLSX from "xlsx";
import { db, auth } from "@/src/lib/firebase";
import { collection, doc, serverTimestamp, getDocs, query, orderBy, limit, getCountFromServer } from "firebase/firestore";
import { signInAnonymously, onAuthStateChanged } from "firebase/auth";
import { logQueryToDatabase } from "@/src/lib/logbook";

import { 
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";

export default function App() {
  const reportRef = useRef<HTMLDivElement>(null);
  const [compounds, setCompounds] = useState<CompoundInput[]>([
    { value: "", type: "SMILES" }
  ]);
  const [view, setView] = useState<'input' | 'loading' | 'results' | 'logbook'>('input');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<PredictionResult | null>(null);
  const [error, setError] = useState<{ message: string; type: string } | null>(null);
  const [user, setUser] = useState<any>(null);
  const [logbookCount, setLogbookCount] = useState<number>(0);

  const fetchLogbookStats = async () => {
    try {
      const snap = await getCountFromServer(collection(db, "query_logs"));
      setLogbookCount(snap.data().count);
    } catch (e) {
      console.warn("Logbook count temporarily unavailable:", e);
    }
  };

  useEffect(() => {
    const init = async () => {
      try {
        initRDKit().catch(console.error);
        await fetchLogbookStats();
      } catch (e) {
        console.error("Initialization Error:", e);
      }
    };

    const unsubscribe = onAuthStateChanged(auth, (u) => {
      setUser(u);
    });

    init();
    return () => unsubscribe();
  }, []);

  const [predictionMethod, setPredictionMethod] = useState<PredictionMethod>("Both");

  const addCompound = () => {
    if (compounds.length < 5) {
      setCompounds([...compounds, { value: "", type: "SMILES" }]);
    }
  };

  const removeCompound = (index: number) => {
    if (compounds.length > 1) {
      const newCompounds = [...compounds];
      newCompounds.splice(index, 1);
      setCompounds(newCompounds);
    }
  };

  const updateCompound = (index: number, field: keyof CompoundInput, value: string) => {
    const newCompounds = [...compounds];
    newCompounds[index] = { ...newCompounds[index], [field]: value, type: "SMILES" };
    setCompounds(newCompounds);
  };

  const handlePredict = async (e: FormEvent) => {
    e.preventDefault();
    const validInputs = compounds.filter(c => c.value.trim() !== "");
    if (validInputs.length === 0) return;

    // Enforce SMILES-only inputs
    const formattedInputs: CompoundInput[] = validInputs.map(c => ({
      value: c.value.trim(),
      type: "SMILES"
    }));

    setView('loading');
    setLoading(true);
    setError(null);
    try {
      // Compute RDKit descriptors for input compounds before passing to AI
      const validInputsWithDescriptors = await Promise.all(formattedInputs.map(async (input) => {
        try {
          const desc = await getMolecularDescriptors(input.value);
          return { ...input, descriptors: desc };
        } catch (e) {
          return input;
        }
      }));

      let switchedView = false;
      const prediction = await predictInteraction(validInputsWithDescriptors, predictionMethod, (partial) => {
        if (partial.degradationImpurities) {
          partial.degradationImpurities = [...partial.degradationImpurities]
            .sort((a, b) => (b.probability || 0) - (a.probability || 0))
            .slice(0, 5);
        }
        if (!switchedView && (
          (partial.chainOfThought && partial.chainOfThought.length > 0) || 
          (partial.compounds && partial.compounds.length > 0)
        )) {
          setView('results');
          setLoading(false);
          switchedView = true;
        }
        setResult(partial as PredictionResult);
      });
      
      if (prediction.degradationImpurities) {
        prediction.degradationImpurities = [...prediction.degradationImpurities]
          .sort((a, b) => (b.probability || 0) - (a.probability || 0))
          .slice(0, 5);
      }

      // Calculate real molecular descriptors (Molecular Weight) using RDKit in parallel
      await Promise.all([
        ...prediction.compounds.map(async (comp) => {
          if (comp.smiles) {
            const descriptors = await getMolecularDescriptors(comp.smiles);
            if (descriptors) {
              comp.molecularDescriptors = descriptors;
            }
          }
        }),
        ...prediction.degradationImpurities.map(async (impurity) => {
          if (impurity.smiles) {
            const descriptors = await getMolecularDescriptors(impurity.smiles);
            if (descriptors) {
              impurity.molecularDescriptors = descriptors;
            }
            if (predictionMethod === "Boltzmann" || predictionMethod === "Both") {
              const strainEnergy = await computeStrainEnergy(impurity.smiles);
              if (strainEnergy !== null) {
                impurity.relativeEnergy = strainEnergy; 
              }
            }
          }
        })
      ]);

      setResult(prediction);
      setView('results');

      // Maintain the CSV logbook of the 100 most recent queries on the database
      const primarySmiles = formattedInputs[0].value;
      const secondarySmilesList = formattedInputs.slice(1).map(c => c.value);
      await logQueryToDatabase(
        primarySmiles,
        secondarySmilesList,
        prediction.degradationImpurities || []
      );

      // Refresh logbook count
      await fetchLogbookStats();

    } catch (err: any) {
      setView('input');
      console.error(err);
      if (err instanceof AnalysisError) {
        setError({ message: err.message, type: err.type });
      } else {
        setError({ 
          message: err.message || "An unexpected error occurred during chemical analysis.", 
          type: "UNKNOWN_ERROR" 
        });
      }
    } finally {
      setLoading(false);
    }
  };

  const downloadExcel = () => {
    if (!result) return;
    
    setLoading(true);
    try {
      const wb = XLSX.utils.book_new();
      const rows: any[][] = [];

      // Title & Header
      rows.push(["INTERACTION REPORT"]);
      rows.push([`Generated on: ${new Date().toLocaleString('en-US', { hour12: false })}`]);
      rows.push([]);

      // 1. Compounds Section
      rows.push(["INPUT COMPOUNDS"]);
      rows.push(["Role", "Name", "SMILES", "MW (g/mol)", "Features", "Interaction Sites"]);
      result.compounds.forEach((c, idx) => {
        rows.push([
          idx === 0 ? "Primary" : "Secondary",
          c.name,
          c.smiles,
          c.molecularDescriptors?.MolWt != null ? c.molecularDescriptors.MolWt.toFixed(2) : "N/A",
          c.features.join(", "),
          c.interactionSites?.join(", ") || "N/A"
        ]);
      });
      rows.push([]);

      // 2. Impurities Section (Top 5 Byproducts)
      const topImpurities = [...(result.degradationImpurities || [])]
        .sort((a, b) => (b.probability || 0) - (a.probability || 0))
        .slice(0, 5);

      if (topImpurities.length > 0) {
        rows.push(["PREDICTED REACTION PRODUCTS & BYPRODUCTS (TOP 5)"]);
        const hasEnergy = topImpurities.some(i => i.relativeEnergy != null);
        const hasBoth = topImpurities.some(i => i.probabilityHeuristic != null);
        
        const header = ["IUPAC Name", "SMILES", "MW (g/mol)", "Main Probability (%)"];
        if (hasBoth) {
          header.push("Heuristic (%)", "Boltzmann (%)");
        }
        if (hasEnergy) header.push("Relative Energy (kcal/mol)");
        header.push("Origin", "Condition", "Source", "Description");
        rows.push(header);

        topImpurities.forEach(i => {
            const row = [
              i.iupacName,
              i.smiles,
              i.molecularDescriptors?.MolWt != null ? i.molecularDescriptors.MolWt.toFixed(2) : "N/A",
              i.probability != null ? (i.probability * 100).toFixed(1) : "N/A"
            ];
            if (hasBoth) {
              row.push(
                i.probabilityHeuristic != null ? (i.probabilityHeuristic * 100).toFixed(1) : "N/A",
                i.probabilityBoltzmann != null ? (i.probabilityBoltzmann * 100).toFixed(1) : "N/A"
              );
            }
            if (hasEnergy) row.push(i.relativeEnergy != null ? i.relativeEnergy.toFixed(2) : "N/A");
            row.push(
              i.origin,
              i.condition,
              i.source,
              i.structureDescription
            );
            rows.push(row);
          });
        rows.push([]);
      }

      // 3. Mechanism Section
      rows.push(["INTERACTION MECHANISM"]);
      rows.push([result.mechanism]);

      const ws = XLSX.utils.aoa_to_sheet(rows);

      // Basic column width adjustments
      const wscols = [
        { wch: 15 }, // Role
        { wch: 25 }, // Name
        { wch: 40 }, // SMILES
        { wch: 30 }, // Features
        { wch: 30 }, // Interaction Sites
        { wch: 20 }, // Condition/Source
        { wch: 50 }, // Description
      ];
      ws['!cols'] = wscols;

      XLSX.utils.book_append_sheet(wb, ws, "Interaction Report");
      XLSX.writeFile(wb, `Interaction_Report_${new Date().toISOString().split('T')[0]}.xlsx`);
    } catch (err) {
      console.error('Excel Generation Error:', err);
      setError({ 
        message: "Failed to generate Excel report. Please try again.", 
        type: "EXCEL_ERROR" 
      });
    } finally {
      setLoading(false);
    }
  };

  const getInteractionIcon = (type: string) => {
    switch (type) {
      case "Chemical": return <ShieldAlert className="w-5 h-5 text-red-500" />;
      case "Physical": return <AlertTriangle className="w-5 h-5 text-amber-500" />;
      default: return <CheckCircle2 className="w-5 h-5 text-green-500" />;
    }
  };

  return (
    <TooltipProvider>
      <div className="min-h-screen bg-white flex flex-col font-sans">
      {/* Header */}
      <header className="border-b border-[#E2E8F0] bg-white sticky top-0 z-20">
        <div className="max-w-[1120px] mx-auto px-4 sm:px-6 h-16 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <h1 className="font-serif text-2xl sm:text-3xl font-extrabold text-[#0F172A] tracking-tight">Interaction</h1>
            {view === 'results' && (
              <span className="text-xs font-semibold px-2 py-0.5 rounded bg-slate-100 text-slate-600">
                Report
              </span>
            )}
          </div>

          <div className="flex items-center gap-2">
            <button
              onClick={() => {
                if (view === 'logbook') {
                  setView(result ? 'results' : 'input');
                }
              }}
              className={`px-3 py-1.5 text-xs font-semibold rounded-lg transition-colors ${
                view !== 'logbook'
                  ? 'bg-[#EEF2FF] text-[#4F46E5]'
                  : 'text-[#64748B] hover:text-[#0F172A] hover:bg-slate-100'
              }`}
            >
              Reaction Modeling
            </button>
            <button
              onClick={() => setView('logbook')}
              className={`inline-flex items-center gap-1.5 px-3 py-1.5 text-xs font-semibold rounded-lg transition-colors ${
                view === 'logbook'
                  ? 'bg-[#EEF2FF] text-[#4F46E5]'
                  : 'text-[#64748B] hover:text-[#0F172A] hover:bg-slate-100'
              }`}
            >
              <FileSpreadsheet className="w-3.5 h-3.5" />
              CSV Logbook
              <span className="text-[10px] px-1.5 py-0.2 rounded-full bg-slate-200 text-slate-700 font-mono">
                {logbookCount}
              </span>
            </button>
            {view === 'results' && (
              <button
                onClick={() => { setView('input'); setResult(null); }}
                className="ml-2 text-xs font-medium text-slate-500 hover:text-slate-900 transition-colors"
              >
                &larr; New Setup
              </button>
            )}
          </div>
        </div>
      </header>

      <main className="flex-1 max-w-[1120px] mx-auto px-4 sm:px-6 py-6 w-full relative">
        {/* CSV Logbook View */}
        {view === 'logbook' && (
          <CsvLogbook 
            onNewReactionClick={() => {
              setView('input');
              setResult(null);
            }} 
          />
        )}

        {/* Loading View Matching Streamlit */}
        {view === 'loading' && (
          <div className="py-20 px-4 text-center max-w-xl mx-auto">
            <div className="inline-block w-14 h-14 border-4 border-[#EEF2FF] border-t-[#4F46E5] rounded-full animate-spin mb-6"></div>
            <h2 className="font-serif text-2xl sm:text-3xl font-bold text-[#0F172A] mb-3">
              Computing Reaction Products & Transformation Pathways...
            </h2>
            <p className="text-[#64748B] text-sm leading-relaxed">
              Analyzing electrophilic and nucleophilic reactive centers, evaluating transition state kinetic activation barriers, and calculating Boltzmann thermodynamic free energies (ΔG).
            </p>
          </div>
        )}
        
        {/* Input View Matching Streamlit */}
        {view === 'input' && (
          <div className="w-full space-y-6">
            <div className="bg-white border border-[#E2E8F0] rounded-2xl p-6 sm:p-8 shadow-[0_1px_3px_rgba(15,23,42,0.03)]">
              <div className="mb-6">
                <h2 className="font-serif text-2xl font-bold text-[#0F172A] mb-1">Reaction Mixture Setup</h2>
                <p className="text-sm text-[#64748B]">
                  Enter Simplified Molecular Input Line Entry System (SMILES) strings for the primary compound and optional secondary co-reactants.
                </p>
              </div>

              {error && (
                <div className="mb-6">
                  <Alert variant="destructive" className="border-red-200 bg-red-50">
                    <div className="flex gap-3">
                      <div className="mt-0.5">
                        {(error.type === "QUOTA_EXCEEDED" || error.type === "MODEL_OVERLOADED") && <Clock className="h-5 w-5 text-red-600 animate-pulse" />}
                        {error.type === "SAFETY_TRIGGERED" && <ShieldAlert className="h-5 w-5 text-red-600" />}
                        {error.type === "CONNECTION_ERROR" && <WifiOff className="h-5 w-5 text-red-600" />}
                        {(error.type === "INVALID_SMILES" || error.type === "INVALID_JSON") && <AlertCircle className="h-5 w-5 text-red-600" />}
                        {(error.type === "CONFIG_ERROR" || error.type === "PERMISSION_DENIED") && <Lock className="h-5 w-5 text-red-600" />}
                        {error.type === "UNKNOWN_ERROR" && <AlertTriangle className="h-5 w-5 text-red-600" />}
                      </div>
                      <div className="space-y-1">
                        <AlertTitle className="text-red-800 font-bold">
                          {error.type === "QUOTA_EXCEEDED" ? "API Quota Exceeded" : 
                           error.type === "MODEL_OVERLOADED" ? "Model Temporarily Overloaded" :
                           error.type === "SAFETY_TRIGGERED" ? "Safety Filter Logic Engaged" :
                           error.type === "INVALID_SMILES" ? "Chemical Structure Error" :
                           error.type === "CONNECTION_ERROR" ? "Network Communication Failure" :
                           error.type === "CONFIG_ERROR" ? "Configuration Credential Error" :
                           error.type === "PERMISSION_DENIED" ? "API Access Permission Denied" :
                           error.type === "INVALID_JSON" ? "Structure Interpretation Failure" :
                           "Analytical Processing Error"}
                        </AlertTitle>
                        <AlertDescription className="text-red-700">
                          {error.message}
                        </AlertDescription>
                        {error.type === "CONFIG_ERROR" && (
                          <div className="mt-3 p-3.5 bg-red-100/80 border border-red-200 rounded-lg text-xs text-red-900 space-y-1.5">
                            <div className="font-bold">How to resolve this in Google AI Studio:</div>
                            <ol className="list-decimal list-inside space-y-1 text-red-800">
                              <li>Obtain a Gemini API key at <a href="https://aistudio.google.com/apikey" target="_blank" rel="noreferrer" className="underline font-semibold hover:text-red-950">aistudio.google.com/apikey</a>.</li>
                              <li>In the AI Studio menu, open <strong>Settings</strong> &gt; <strong>Secrets</strong> (or Environment Variables).</li>
                              <li>Update <strong>GEMINI_API_KEY</strong> with your genuine Gemini API key (note that Firebase Web API keys cannot be used for Generative AI).</li>
                            </ol>
                          </div>
                        )}
                        {(error.type === "QUOTA_EXCEEDED" || error.type === "MODEL_OVERLOADED" || error.type === "CONNECTION_ERROR" || error.type === "PERMISSION_DENIED" || error.type === "UNKNOWN_ERROR" || error.type === "CONFIG_ERROR") && (
                          <Button 
                            variant="outline" 
                            size="sm" 
                            onClick={(e) => handlePredict(e as any)}
                            className="mt-3 border-red-200 text-red-700 hover:bg-red-100"
                          >
                            <RefreshCw className="mr-2 h-3 w-3" />
                            Retry Analytical Operation
                          </Button>
                        )}
                      </div>
                    </div>
                  </Alert>
                </div>
              )}

              <form onSubmit={handlePredict} autoComplete="off">
                {/* Section: Primary Compound (SMILES Only) */}
                <div className="mb-6">
                  <div className="flex items-center justify-between mb-2">
                    <div className="flex items-center gap-2 text-sm font-bold text-[#312E81]">
                      <span className="inline-block w-2 h-2 rounded-full bg-[#4F46E5]"></span>
                      Primary Compound (SMILES)
                    </div>
                    <span className="text-[11px] font-medium text-[#64748B]">SMILES required</span>
                  </div>
                  <div>
                    <input 
                      placeholder="e.g. CC(=O)Oc1ccccc1C(=O)O (Aspirin)"
                      value={compounds[0].value}
                      autoComplete="off"
                      spellCheck={false}
                      onChange={(e) => {
                        updateCompound(0, "value", e.target.value);
                        setError(null);
                      }}
                      required
                      className="w-full h-11 px-3.5 font-mono text-sm bg-white border border-[#E2E8F0] rounded-lg text-[#0F172A] placeholder:text-[#94A3B8] placeholder:font-sans focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
                    />
                    <p className="mt-1.5 text-xs text-[#64748B]">
                      SMILES only. E.g. Aspirin: <code className="font-mono bg-slate-100 px-1.5 py-0.5 rounded text-[#0F172A]">CC(=O)Oc1ccccc1C(=O)O</code>
                    </p>
                  </div>
                </div>

                {/* Section: Secondary Compounds */}
                <div className="mb-6">
                  <div className="flex items-center justify-between mb-3">
                    <div className="flex items-center gap-2 text-sm font-bold text-[#334155]">
                      <span className="inline-block w-2 h-2 rounded-full bg-[#94A3B8]"></span>
                      Secondary Compounds (SMILES) (Co-reactants / Additives)
                    </div>
                    <span className="font-mono text-xs text-[#64748B] bg-[#F1F5F9] px-2 py-0.5 rounded">
                      {compounds.slice(1).filter(c => c.value.trim()).length} Added
                    </span>
                  </div>

                  <div className="space-y-3 mb-4">
                    {compounds.slice(1).map((c, idx) => {
                      const actualIndex = idx + 1;
                      return (
                        <div key={`sec-${actualIndex}`} className="flex items-center gap-3">
                          <div className="flex-1">
                            <input
                              placeholder={`e.g. Secondary Compound ${actualIndex} SMILES (e.g. CC(=O)NC1=CC=C(O)C=C1)`}
                              value={c.value}
                              autoComplete="off"
                              spellCheck={false}
                              onChange={(e) => {
                                updateCompound(actualIndex, "value", e.target.value);
                                setError(null);
                              }}
                              className="w-full h-10 px-3.5 font-mono text-sm bg-white border border-[#E2E8F0] rounded-lg text-[#0F172A] placeholder:text-[#94A3B8] placeholder:font-sans focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
                            />
                          </div>
                          <button
                            type="button"
                            onClick={() => removeCompound(actualIndex)}
                            className="px-3 h-10 flex items-center justify-center text-[#94A3B8] hover:text-red-600 hover:bg-red-50 border border-[#E2E8F0] hover:border-red-200 rounded-lg transition-colors text-xs font-medium"
                          >
                            Remove
                          </button>
                        </div>
                      );
                    })}
                  </div>

                  {compounds.length < 5 && (
                    <button
                      type="button"
                      onClick={addCompound}
                      className="inline-flex items-center gap-1.5 px-4 py-2 border border-[#CBD5E1] bg-white text-[#334155] hover:bg-[#F8FAFC] text-xs font-semibold rounded-lg shadow-xs transition-colors"
                    >
                      + Add Secondary Compound (SMILES)
                    </button>
                  )}
                </div>

                {/* Section: Prediction Method */}
                <div className="mb-6">
                  <div className="text-xs font-bold text-[#64748B] uppercase tracking-wider mb-2.5">
                    Prediction Engine & Methodology:
                  </div>
                  <div className="space-y-2">
                    {[
                      { id: "Both", title: "Dual Engine (Heuristic Kinetic Rules + Boltzmann Thermodynamic ΔG)" },
                      { id: "Heuristic", title: "Heuristic (Expert Kinetic Activation & Transition States)" },
                      { id: "Boltzmann", title: "Boltzmann (Thermodynamic Free Energy ΔG Distribution at 298.15K)" }
                    ].map((opt) => (
                      <label
                        key={opt.id}
                        onClick={() => setPredictionMethod(opt.id as PredictionMethod)}
                        className={`flex items-center gap-3 p-3 rounded-xl border cursor-pointer transition-all ${
                          predictionMethod === opt.id
                            ? "bg-[#EEF2FF] border-[#818CF8] text-[#312E81] shadow-xs"
                            : "bg-white border-[#E2E8F0] hover:border-[#CBD5E1] text-[#334155]"
                        }`}
                      >
                        <div className={`w-4 h-4 rounded-full border-2 flex items-center justify-center ${
                          predictionMethod === opt.id ? "border-[#4F46E5] bg-[#4F46E5]" : "border-[#94A3B8] bg-white"
                        }`}>
                          {predictionMethod === opt.id && <div className="w-1.5 h-1.5 bg-white rounded-full" />}
                        </div>
                        <span className="text-xs sm:text-sm font-semibold">{opt.title}</span>
                      </label>
                    ))}
                  </div>
                </div>

                {/* Submit CTA Button */}
                <button
                  type="submit"
                  disabled={loading || compounds.every(c => c.value.trim() === "")}
                  className="w-full py-3 px-6 bg-[#4F46E5] hover:bg-[#4338CA] disabled:opacity-50 text-white font-semibold text-base rounded-lg shadow-sm hover:shadow transition-all flex items-center justify-center gap-2"
                >
                  Predict Chemical Interactions
                </button>
              </form>
            </div>
          </div>
        )}

        {/* View: Results Dashboard Matching Streamlit */}
        {view === 'results' && result && !loading && (() => {
          // Explicitly constrain degradationImpurities to only the top 5 elements ranked by formation probability
          const topImpurities = [...(result.degradationImpurities || [])]
            .sort((a, b) => (b.probability || 0) - (a.probability || 0))
            .slice(0, 5);

          return (
            <div className="w-full space-y-6">
            {/* 1. Input Chemical Data Card */}
            <section className="bg-white border border-[#E2E8F0] rounded-2xl p-6 shadow-xs">
              <div className="mb-5">
                <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                  Input Chemical Data
                </h3>
                <p className="text-xs sm:text-sm text-[#64748B]">
                  Calculated molecular descriptors, functional group features, and predicted reactive interaction sites.
                </p>
              </div>

              <div className="flex flex-col gap-4">
                {(result.compounds || []).map((comp, idx) => {
                  const mw = comp.molecularDescriptors?.MolWt;
                  const role = idx === 0 ? "Primary Compound" : `Secondary Compound ${idx}`;
                  const roleClass = idx === 0 ? "role-primary" : "role-secondary";

                  return (
                    <div key={`comp-card-${idx}`} className="ap1-comp-card flex-col sm:flex-row w-full">
                      <div className="ap1-comp-mol">
                        <span className="ap1-comp-badge">C{idx + 1}</span>
                        {comp.smiles ? (
                          <ChemicalStructure smiles={comp.smiles} width={180} height={180} />
                        ) : (
                          <div className="w-36 h-36 bg-slate-100 rounded-lg animate-pulse" />
                        )}
                      </div>
                      <div className="ap1-comp-info">
                        <div className="flex items-center mb-1">
                          <span className="ap1-comp-name">{comp.name || "Compound"}</span>
                          <span className={`ap1-comp-role ${roleClass}`}>{role}</span>
                        </div>
                        {comp.smiles && (
                          <div className="ap1-smiles-box" title={comp.smiles}>
                            {comp.smiles}
                          </div>
                        )}
                        <div className="ap1-tag-group">
                          {mw && <span className="ap1-pill mw">MW: {mw.toFixed(2)} g/mol</span>}
                          {(comp.features || []).map((f, fi) => (
                            <span key={fi} className="ap1-pill">{f}</span>
                          ))}
                        </div>
                        {comp.interactionSites && comp.interactionSites.length > 0 && (
                          <div className="mt-3">
                            <div className="text-[11px] font-bold text-[#2563EB] uppercase tracking-wider mb-1">
                              Reactive Interaction Centers:
                            </div>
                            <div className="ap1-tag-group">
                              {comp.interactionSites.map((site, si) => (
                                <span key={si} className="ap1-pill site">{site}</span>
                              ))}
                            </div>
                          </div>
                        )}
                      </div>
                    </div>
                  );
                })}
              </div>
            </section>

            {/* Functional Group Reactivity Analysis Card */}
            {(() => {
              const primaryComp = (result.compounds || [])[0];
              const fgList: FunctionalGroupReactivity[] = result.functionalGroupAnalysis || 
                (primaryComp?.smiles ? detectFunctionalGroupsDetailed(primaryComp.smiles).functionalGroups : []);
              
              if (fgList.length === 0) return null;

              const getVulnStyle = (vuln: string) => {
                switch (vuln) {
                  case "Critical": return "bg-rose-50 text-rose-700 border-rose-200";
                  case "High": return "bg-amber-50 text-amber-700 border-amber-200";
                  case "Moderate": return "bg-yellow-50 text-yellow-800 border-yellow-200";
                  case "Low": return "bg-blue-50 text-blue-700 border-blue-200";
                  default: return "bg-slate-50 text-slate-600 border-slate-200";
                }
              };

              return (
                <section className="bg-white border border-[#E2E8F0] rounded-2xl p-6 shadow-xs">
                  <div className="mb-5">
                    <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                      Functional Group Reactivity Analysis
                    </h3>
                    <p className="text-xs sm:text-sm text-[#64748B]">
                      Systematic evaluation of identified functional groups and their mechanistic reactivity with acidic, basic, hydrolysis, photolytic, thermal, oxidative conditions, and co-reactant functional groups.
                    </p>
                  </div>

                  <div className="flex flex-col gap-5">
                    {fgList.map((fg, fgIdx) => (
                      <div key={`fg-${fgIdx}`} className="border border-[#E2E8F0] rounded-xl p-4 bg-[#FAFAFA] flex flex-col gap-3">
                        <div className="flex flex-wrap items-center justify-between gap-2 pb-2 border-b border-[#E2E8F0]">
                          <div className="flex items-center gap-2">
                            <span className="font-semibold text-sm text-[#0F172A]">{fg.groupName}</span>
                            <span className="text-[11px] font-medium px-2 py-0.5 rounded bg-slate-200 text-slate-700">
                              {fg.category}
                            </span>
                            <span className="font-mono text-xs px-2 py-0.5 rounded bg-indigo-50 text-indigo-700 border border-indigo-100">
                              {fg.smilesFragment}
                            </span>
                          </div>
                          <div className="text-xs text-[#64748B]">
                            <span className="font-medium text-[#475569]">Reactive Center:</span> {fg.reactiveSite}
                          </div>
                        </div>

                        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
                          {/* Acidic */}
                          <div className="bg-white border border-[#E2E8F0] rounded-lg p-3">
                            <div className="flex items-center justify-between mb-1.5">
                              <span className="text-xs font-bold text-[#0F172A]">Acidic Stress</span>
                              <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.acidic.vulnerability)}`}>
                                {fg.acidic.vulnerability}
                              </span>
                            </div>
                            <p className="text-[11px] text-[#475569] leading-relaxed">{fg.acidic.mechanism}</p>
                          </div>

                          {/* Basic */}
                          <div className="bg-white border border-[#E2E8F0] rounded-lg p-3">
                            <div className="flex items-center justify-between mb-1.5">
                              <span className="text-xs font-bold text-[#0F172A]">Basic Stress</span>
                              <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.basic.vulnerability)}`}>
                                {fg.basic.vulnerability}
                              </span>
                            </div>
                            <p className="text-[11px] text-[#475569] leading-relaxed">{fg.basic.mechanism}</p>
                          </div>

                          {/* Hydrolysis */}
                          <div className="bg-white border border-[#E2E8F0] rounded-lg p-3">
                            <div className="flex items-center justify-between mb-1.5">
                              <span className="text-xs font-bold text-[#0F172A]">Hydrolysis (Aqueous)</span>
                              <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.hydrolysis.vulnerability)}`}>
                                {fg.hydrolysis.vulnerability}
                              </span>
                            </div>
                            <p className="text-[11px] text-[#475569] leading-relaxed">{fg.hydrolysis.mechanism}</p>
                          </div>

                          {/* Photolytic */}
                          <div className="bg-white border border-[#E2E8F0] rounded-lg p-3">
                            <div className="flex items-center justify-between mb-1.5">
                              <span className="text-xs font-bold text-[#0F172A]">Photolytic Stress</span>
                              <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.photolytic.vulnerability)}`}>
                                {fg.photolytic.vulnerability}
                              </span>
                            </div>
                            <p className="text-[11px] text-[#475569] leading-relaxed">{fg.photolytic.mechanism}</p>
                          </div>

                          {/* Thermal */}
                          <div className="bg-white border border-[#E2E8F0] rounded-lg p-3">
                            <div className="flex items-center justify-between mb-1.5">
                              <span className="text-xs font-bold text-[#0F172A]">Thermal Stress</span>
                              <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.thermal.vulnerability)}`}>
                                {fg.thermal.vulnerability}
                              </span>
                            </div>
                            <p className="text-[11px] text-[#475569] leading-relaxed">{fg.thermal.mechanism}</p>
                          </div>

                          {/* Oxidative */}
                          <div className="bg-white border border-[#E2E8F0] rounded-lg p-3">
                            <div className="flex items-center justify-between mb-1.5">
                              <span className="text-xs font-bold text-[#0F172A]">Oxidative Stress</span>
                              <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.oxidative.vulnerability)}`}>
                                {fg.oxidative.vulnerability}
                              </span>
                            </div>
                            <p className="text-[11px] text-[#475569] leading-relaxed">{fg.oxidative.mechanism}</p>
                          </div>

                          {/* Secondary Interaction (if present) */}
                          {fg.secondaryInteraction && (
                            <div className="bg-white border border-indigo-200 rounded-lg p-3 md:col-span-2 lg:col-span-3">
                              <div className="flex items-center justify-between mb-1.5">
                                <span className="text-xs font-bold text-indigo-900">
                                  Secondary Compound Cross-Reactivity {fg.secondaryInteraction.partnerGroup ? `(Target: ${fg.secondaryInteraction.partnerGroup})` : ""}
                                </span>
                                <span className={`text-[11px] font-semibold px-2 py-0.5 rounded border ${getVulnStyle(fg.secondaryInteraction.vulnerability)}`}>
                                  {fg.secondaryInteraction.vulnerability}
                                </span>
                              </div>
                              <p className="text-[11px] text-[#475569] leading-relaxed">{fg.secondaryInteraction.mechanism}</p>
                            </div>
                          )}
                        </div>
                      </div>
                    ))}
                  </div>
                </section>
              );
            })()}

            {/* 2. Interaction Potential & Reactive Centers Heatmap (Seaborn) */}
            <InteractionHeatmap compounds={result.compounds || []} />

            {/* 3. Mechanistic Framework Evaluation Card */}
            <section className="bg-white border border-[#E2E8F0] rounded-2xl p-6 shadow-xs">
              <div className="mb-5">
                <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                  Mechanistic Framework Evaluation
                </h3>
                <p className="text-xs sm:text-sm text-[#64748B]">
                  Comprehensive kinetic pathways, microenvironmental influences, and thermodynamic justification.
                </p>
              </div>

              <div className="flex flex-col gap-4">
                <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl p-5 text-sm text-[#334155] leading-relaxed whitespace-pre-wrap font-sans">
                  {result.chainOfThought || "No detailed reasoning chain provided."}
                </div>

                <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl p-4 text-xs text-[#475569] leading-relaxed">
                  <strong className="text-[#0F172A] block mb-1 text-sm font-semibold">
                    Chemical Reaction & Byproduct Analysis:
                  </strong>
                  Products identified with high formation probability or favorable exergonic free energy (ΔG &lt; 0 kcal/mol) represent dominant reaction pathways. In experimental validation, these byproducts should be verified using analytical separation techniques (HPLC, LC-MS, GC-MS, or NMR).
                </div>
              </div>
            </section>

            {/* 4. Degradation Products and Details Card */}
            <section className="bg-white border border-[#E2E8F0] rounded-2xl p-6 shadow-xs">
              <div className="mb-5">
                <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                  Degradation Products and Details
                </h3>
                <p className="text-xs sm:text-sm text-[#64748B]">
                  Ranked strictly by formation probability and thermodynamic stability (Top 5 maximum).
                </p>
              </div>

              {topImpurities.length === 0 ? (
                <div className="p-8 text-center bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl text-sm text-[#64748B]">
                  No significant byproducts detected under standard conditions.
                </div>
              ) : (
                <div className="flex flex-col gap-4">
                  {topImpurities.map((imp, idx) => {
                    const prob = (imp.probability || 0) * 100;
                    const cond = imp.condition || "Direct Degradation";
                    let condClass = "cond-hydro";
                    const condLower = cond.toLowerCase();
                    if (condLower.includes("oxid")) condClass = "cond-oxid";
                    else if (condLower.includes("therm")) condClass = "cond-therm";
                    else if (condLower.includes("photo")) condClass = "cond-photo";
                    else if (condLower.includes("react") || condLower.includes("incomp")) condClass = "cond-react";

                    return (
                      <div key={`prod-${idx}`} className="ap1-imp-card flex-col sm:flex-row w-full">
                        <div className="ap1-imp-svg">
                          <div className="absolute top-3 left-3 bg-[#F1F5F9] text-[#475569] text-xs font-extrabold px-2 py-0.5 rounded">
                            #{idx + 1}
                          </div>
                          {imp.smiles ? (
                            <ChemicalStructure smiles={imp.smiles} width={220} height={220} />
                          ) : (
                            <div className="w-40 h-40 bg-slate-100 rounded-lg animate-pulse" />
                          )}
                        </div>
                        <div className="ap1-imp-body">
                          <div className="ap1-imp-header">
                            <div>
                              <div className="ap1-imp-title">{imp.iupacName || "Product"}</div>
                              <div className="mt-1.5 flex gap-2">
                                {imp.molecularDescriptors?.MolWt && (
                                  <span className="ap1-pill mw">MW: {imp.molecularDescriptors.MolWt.toFixed(2)} g/mol</span>
                                )}
                              </div>
                            </div>
                            <div className="text-right">
                              <div className="ap1-imp-prob-val">{prob.toFixed(1)}%</div>
                              {imp.probabilityHeuristic != null && imp.probabilityBoltzmann != null && (
                                <div className="ap1-imp-prob-sub">
                                  Heuristic: {(imp.probabilityHeuristic * 100).toFixed(1)}% | Boltzmann: {(imp.probabilityBoltzmann * 100).toFixed(1)}%
                                </div>
                              )}
                              {imp.relativeEnergy != null && (
                                <div className="text-xs font-mono text-[#64748B] text-right mt-1">
                                  ΔG: {imp.relativeEnergy.toFixed(2)} kcal/mol
                                </div>
                              )}
                            </div>
                          </div>

                          <div className="ap1-prob-bar-bg">
                            <div className="ap1-prob-bar-fill" style={{ width: `${Math.min(Math.max(prob, 5), 100)}%` }}></div>
                          </div>

                          {imp.structureDescription && (
                            <div className="ap1-imp-desc">{imp.structureDescription}</div>
                          )}

                          {imp.mechanismExplanation && (
                            <div className="ap1-mech-box">
                              <div className="ap1-mech-title">
                                <span>Chemical Mechanism:</span>
                              </div>
                              <div>{imp.mechanismExplanation}</div>
                            </div>
                          )}

                          <div className="flex flex-wrap gap-2 items-center">
                            <span className={`ap1-badge-cond ${condClass}`}>{cond}</span>
                            <span className="ap1-pill font-semibold text-[#1D4ED8] bg-[#EFF6FF] border-[#DBEAFE]">
                              Origin: {imp.origin || "Parent Molecule"}
                            </span>
                            {imp.smiles && (
                              <span className="ap1-pill font-mono text-[11px] text-[#64748B] truncate max-w-xs" title={imp.smiles}>
                                {imp.smiles}
                              </span>
                            )}
                          </div>
                        </div>
                      </div>
                    );
                  })}
                </div>
              )}
            </section>

            {/* 5. Disclaimer Card */}
            <div className="border border-[#E2E8F0] rounded-xl bg-[#FAFAFA] p-4 text-xs text-[#64748B] leading-relaxed">
              <p className="italic m-0">
                Disclaimer: INTERACTION is an AI-assisted computational chemistry modeling tool designed for reaction pathway exploration and byproduct screening. Predictions should be verified by experimental analytical assays (HPLC, LC-MS, NMR).
              </p>
            </div>
          </div>
          );
        })()}
      </main>

      {/* Footer Matching Streamlit */}
      <footer className="border-t border-[#E2E8F0] py-6 bg-white mt-auto">
        <div className="max-w-[1120px] mx-auto px-4 sm:px-6 flex flex-col sm:flex-row justify-between items-center gap-3 text-xs text-[#64748B]">
          <p>© 2026 INTERACTION Chemical Informatics. All rights reserved.</p>
          <div className="flex items-center gap-4 text-[11px]">
            <span>Thermodynamic & Kinetic Modeling</span>
            <span>•</span>
            <span>Computational Chemoinformatics</span>
          </div>
        </div>
      </footer>
    </div>
    </TooltipProvider>
  );
}
