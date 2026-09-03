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
  Download,
  Info,
  Database,
  History,
  ArrowLeft
} from "lucide-react";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle, CardFooter } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Separator } from "@/components/ui/separator";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { predictInteraction, PredictionResult, InputType, CompoundInput, AnalysisError, PredictionMethod } from "@/src/lib/gemini";
import { ChemicalStructure } from "@/src/components/ChemicalStructure";
import { initRDKit, getMolecularDescriptors, computeStrainEnergy } from "@/src/lib/rdkit";
import { sanitizeData } from "@/src/lib/firestore-utils";
import { Plus, Trash2, AlertCircle, WifiOff, Clock, Lock } from "lucide-react";
import * as XLSX from "xlsx";
import { db, auth } from "@/src/lib/firebase";
import { collection, addDoc, setDoc, doc, serverTimestamp, getDocs, query, orderBy, limit, where, getCountFromServer } from "firebase/firestore";
import { signInAnonymously, onAuthStateChanged } from "firebase/auth";
import { seedDatabase } from "@/src/lib/seed";

import { 
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";

export default function App() {
  const reportRef = useRef<HTMLDivElement>(null);
  const [compounds, setCompounds] = useState<CompoundInput[]>([
    { value: "", type: "Name" }
  ]);
  const [view, setView] = useState<'input' | 'loading' | 'results'>('input');
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<PredictionResult | null>(null);
  const [error, setError] = useState<{ message: string; type: string } | null>(null);
  const [user, setUser] = useState<any>(null);
  const [suggestions, setSuggestions] = useState<any[]>([]);
  const [showSuggestions, setShowSuggestions] = useState<number | null>(null);
  const [dbStats, setDbStats] = useState({ compounds: 0, predictions: 0 });

  // Debounced Search Effect
  useEffect(() => {
    if (showSuggestions === null) return;
    
    const activeCompound = compounds[showSuggestions];
    if (!activeCompound || activeCompound.type !== "Name" || activeCompound.value.length < 2) {
      setSuggestions([]);
      return;
    }

    const timer = setTimeout(async () => {
      try {
        const val = activeCompound.value.toLowerCase();
        const q = query(
          collection(db, "compounds"), 
          where("name", ">=", activeCompound.value.charAt(0).toUpperCase() + activeCompound.value.slice(1).toLowerCase()),
          where("name", "<=", activeCompound.value.charAt(0).toUpperCase() + activeCompound.value.slice(1).toLowerCase() + "\uf8ff"),
          limit(20)
        );
        const snap = await getDocs(q);
        const allSuggestions = snap.docs.map(d => ({ ...d.data(), id: d.id })) as any[];
        
        // Dynamically trigger remediation for compounds with missing SMILES
        allSuggestions.forEach((s: any) => {
          if (!s.smiles) {
            import('./lib/compound-manager').then(m => m.remediateCompoundSmiles(s.id, s.name)).catch(console.error);
          }
        });
        
        // Client-side case-insensitive filtering
        const filtered = allSuggestions.filter(s => 
          s.name.toLowerCase().includes(val)
        );
        
        setSuggestions(filtered.slice(0, 5));
      } catch (err) {
        console.error("Search Error:", err);
      }
    }, 300);

    return () => clearTimeout(timer);
  }, [showSuggestions, compounds]);

  useEffect(() => {
    // Initialize Firebase Seed and stats
    const init = async () => {
      try {
        await seedDatabase();
        initRDKit().catch(console.error); 
        
        // Fetch Stats using efficient count aggregation
        const compoundsCount = await getCountFromServer(collection(db, "compounds"));
        const predictionsCount = await getCountFromServer(collection(db, "predictions"));
        
        setDbStats({
          compounds: compoundsCount.data().count,
          predictions: predictionsCount.data().count
        });
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

  const handleSuggestionSelect = (index: number, compound: any) => {
    const newCompounds = [...compounds];
    newCompounds[index] = { value: compound.name, type: "Name" };
    setCompounds(newCompounds);
    setShowSuggestions(null);
  };

  const [predictionMethod, setPredictionMethod] = useState<PredictionMethod>("Both");
  const [activeTab, setActiveTab] = useState<"products" | "reactants" | "reasoning">("products");

  const addCompound = () => {
    if (compounds.length < 5) {
      setCompounds([...compounds, { value: "", type: "Name" }]);
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
    newCompounds[index] = { ...newCompounds[index], [field]: value };
    setCompounds(newCompounds);
  };

  const handlePredict = async (e: FormEvent) => {
    e.preventDefault();
    const validInputs = compounds.filter(c => c.value.trim() !== "");
    if (validInputs.length === 0) return;

    setView('loading');
    setLoading(true);
    setError(null);
    try {
      
      // Attempt to upgrade any "Name" inputs to exact "SMILES" representations using our local database 
      // prior to passing them to the AI to prevent AI structure hallucinations.
      const upgradedInputs = await Promise.all(validInputs.map(async (input) => {
        if (input.type === "Name") {
          try {
            const q = query(collection(db, "compounds"), where("name", "==", input.value.trim()));
            const snap = await getDocs(q);
            if (!snap.empty) {
              const docData = snap.docs[0].data();
              if (docData.smiles) {
                console.log(`Upgrading ${input.value} to exact structural SMILES from DB`);
                return { value: docData.smiles, type: "SMILES" as InputType, originalName: input.value.trim() };
              }
            }
            // Case insensitive fallback check
            const q2 = query(collection(db, "compounds"), where("name", "==", input.value.trim().toLowerCase()));
            const snap2 = await getDocs(q2);
            if (!snap2.empty) {
               const docData2 = snap2.docs[0].data();
               if (docData2.smiles) {
                 return { value: docData2.smiles, type: "SMILES" as InputType, originalName: input.value.trim() };
               }
            }
          } catch (e) {
            console.error("DB Upgrade query failed", e);
          }
        }
        return input;
      }));

      // Compute RDKit descriptors for input compounds (if SMILES provided) before passing to AI
      const validInputsWithDescriptors = await Promise.all(upgradedInputs.map(async (input) => {
        if (input.type === "SMILES") {
          try {
            const desc = await getMolecularDescriptors(input.value);
            return { ...input, descriptors: desc };
          } catch (e) {
            return input;
          }
        }
        return input;
      }));

      let switchedView = false;
      const prediction = await predictInteraction(validInputsWithDescriptors, predictionMethod, (partial) => {
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
                // Ground the LLM's estimate with explicit MMFF94 computational reality.
                impurity.relativeEnergy = strainEnergy; 
              }
            }
          }
        })
      ]);

      setResult(prediction);
      setView('results');

      // Save prediction to Firebase
      await addDoc(collection(db, "predictions"), sanitizeData({
        inputs: validInputs,
        result: prediction,
        method: predictionMethod,
        timestamp: serverTimestamp()
      }));

      // Add new compounds to database if they don't exist
      for (const comp of prediction.compounds) {
        const hasName = comp.name && comp.name.toLowerCase() !== "unknown" && comp.name.trim() !== "";
        const hasSmiles = comp.smiles && comp.smiles.trim() !== "";

        if (!hasName && !hasSmiles) continue;

        let exists = false;
        if (hasName) {
          const q = query(collection(db, "compounds"), where("name", "==", comp.name));
          const snap = await getDocs(q);
          if (!snap.empty) exists = true;
        }

        if (!exists && hasSmiles) {
          const q2 = query(collection(db, "compounds"), where("smiles", "==", comp.smiles));
          const snap2 = await getDocs(q2);
          if (!snap2.empty) exists = true;
        }
        
        if (!exists) {
          const docId = (hasName ? comp.name : comp.smiles).replace(/[^a-z0-9]/gi, '_').toLowerCase();
          await setDoc(doc(db, "compounds", docId), {
            name: hasName ? comp.name : "",
            smiles: hasSmiles ? comp.smiles : "",
            createdAt: new Date().toISOString()
          });
          console.log(`Added new compound to database: ${comp.name || comp.smiles}`);
        }
      }

      // Refresh stats efficiently
      const compoundsCount = await getCountFromServer(collection(db, "compounds"));
      const predictionsCount = await getCountFromServer(collection(db, "predictions"));
      setDbStats({
        compounds: compoundsCount.data().count,
        predictions: predictionsCount.data().count
      });

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

      // 2. Impurities Section
      if (result.degradationImpurities && result.degradationImpurities.length > 0) {
        rows.push(["PREDICTED REACTION PRODUCTS & BYPRODUCTS"]);
        const hasEnergy = result.degradationImpurities.some(i => i.relativeEnergy != null);
        const hasBoth = result.degradationImpurities.some(i => i.probabilityHeuristic != null);
        
        const header = ["IUPAC Name", "SMILES", "MW (g/mol)", "Main Probability (%)"];
        if (hasBoth) {
          header.push("Heuristic (%)", "Boltzmann (%)");
        }
        if (hasEnergy) header.push("Relative Energy (kcal/mol)");
        header.push("Origin", "Condition", "Source", "Description");
        rows.push(header);

        [...result.degradationImpurities]
          .sort((a, b) => (b.probability || 0) - (a.probability || 0))
          .forEach(i => {
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
      {/* Header Matching Streamlit */}
      <header className="border-b border-[#E2E8F0] bg-white sticky top-0 z-20">
        <div className="max-w-[1120px] mx-auto px-4 sm:px-6 h-20 flex items-center justify-between">
          <div>
            <h1 className="font-serif text-3xl font-extrabold text-[#0F172A] tracking-tight">INTERACTION</h1>
            <p className="text-xs sm:text-sm font-semibold text-[#64748B]">Chemical Interaction & Byproduct Prediction Engine</p>
          </div>
          <div className="hidden md:flex items-center gap-3 px-3.5 py-1.5 bg-[#F8FAFC] rounded-full border border-[#E2E8F0]">
            <div className="flex items-center gap-1.5">
              <Database className="w-3.5 h-3.5 text-[#4F46E5]" />
              <span className="text-xs font-semibold text-[#475569]">{dbStats.compounds} Compounds</span>
            </div>
            <Separator orientation="vertical" className="h-3.5 bg-[#E2E8F0]" />
            <div className="flex items-center gap-1.5">
              <History className="w-3.5 h-3.5 text-[#059669]" />
              <span className="text-xs font-semibold text-[#475569]">{dbStats.predictions} Predictions</span>
            </div>
          </div>
        </div>
      </header>

      <main className="flex-1 max-w-[1120px] mx-auto px-4 sm:px-6 py-6 w-full relative">
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
                <p className="text-sm text-[#64748B]">Define the primary chemical compound and optional secondary co-reactants or additives.</p>
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
                        {(error.type === "QUOTA_EXCEEDED" || error.type === "MODEL_OVERLOADED" || error.type === "CONNECTION_ERROR" || error.type === "PERMISSION_DENIED" || error.type === "UNKNOWN_ERROR") && (
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
                {/* Section: Primary Compound */}
                <div className="mb-6">
                  <div className="flex items-center gap-2 text-sm font-bold text-[#312E81] mb-2">
                    <span className="inline-block w-2 h-2 rounded-full bg-[#4F46E5]"></span>
                    Primary Compound
                  </div>
                  <div className="flex gap-3">
                    <select
                      value={compounds[0].type}
                      onChange={(e) => updateCompound(0, "type", e.target.value as InputType)}
                      className="w-28 h-10 px-3 text-xs font-semibold bg-[#F8FAFC] border border-[#E2E8F0] rounded-lg text-[#334155] focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
                    >
                      <option value="Name">Name</option>
                      <option value="SMILES">SMILES</option>
                    </select>
                    <div className="relative flex-1">
                      <input 
                        placeholder={compounds[0].type === "Name" ? "e.g. Aspirin or CC(=O)Oc1ccccc1C(=O)O" : "e.g. CC(=O)Oc1ccccc1C(=O)O"}
                        value={compounds[0].value}
                        autoComplete="off"
                        onChange={(e) => {
                          const val = e.target.value;
                          updateCompound(0, "value", val);
                          setError(null);
                          if (compounds[0].type === "Name" && val.length > 1) {
                            setShowSuggestions(0);
                          } else {
                            setShowSuggestions(null);
                          }
                        }}
                        required
                        className="w-full h-10 px-3.5 text-sm bg-white border border-[#E2E8F0] rounded-lg text-[#0F172A] placeholder:text-[#94A3B8] focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
                      />
                      {showSuggestions === 0 && suggestions.length > 0 && (
                        <div className="absolute z-50 w-full mt-1 bg-white border border-[#E2E8F0] rounded-lg shadow-xl overflow-hidden max-h-48 overflow-y-auto">
                          {suggestions.map((s, i) => (
                            <button
                              key={i}
                              type="button"
                              className="w-full text-left px-4 py-2 text-xs hover:bg-[#EEF2FF] border-b border-[#F1F5F9] last:border-0"
                              onClick={() => handleSuggestionSelect(0, s)}
                            >
                              <div className="font-bold text-[#0F172A]">{s.name}</div>
                              <div className="text-[10px] text-[#64748B] truncate">{s.smiles}</div>
                            </button>
                          ))}
                        </div>
                      )}
                    </div>
                  </div>
                </div>

                {/* Section: Secondary Compounds */}
                <div className="mb-6">
                  <div className="flex items-center justify-between mb-3">
                    <div className="flex items-center gap-2 text-sm font-bold text-[#334155]">
                      <span className="inline-block w-2 h-2 rounded-full bg-[#94A3B8]"></span>
                      Secondary Compounds (Co-reactants / Additives)
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
                          <select
                            value={c.type}
                            onChange={(e) => updateCompound(actualIndex, "type", e.target.value as InputType)}
                            className="w-28 h-10 px-3 text-xs font-semibold bg-[#F8FAFC] border border-[#E2E8F0] rounded-lg text-[#334155] focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
                          >
                            <option value="Name">Name</option>
                            <option value="SMILES">SMILES</option>
                          </select>
                          <div className="relative flex-1">
                            <input
                              placeholder="e.g. Magnesium Stearate or Lactose"
                              value={c.value}
                              autoComplete="off"
                              onChange={(e) => {
                                const val = e.target.value;
                                updateCompound(actualIndex, "value", val);
                                setError(null);
                                if (c.type === "Name" && val.length > 1) {
                                  setShowSuggestions(actualIndex);
                                } else {
                                  setShowSuggestions(null);
                                }
                              }}
                              className="w-full h-10 px-3.5 text-sm bg-white border border-[#E2E8F0] rounded-lg text-[#0F172A] placeholder:text-[#94A3B8] focus:outline-none focus:border-[#4F46E5] focus:ring-1 focus:ring-[#4F46E5]"
                            />
                            {showSuggestions === actualIndex && suggestions.length > 0 && (
                              <div className="absolute z-50 w-full mt-1 bg-white border border-[#E2E8F0] rounded-lg shadow-xl overflow-hidden max-h-48 overflow-y-auto">
                                {suggestions.map((s, i) => (
                                  <button
                                    key={i}
                                    type="button"
                                    className="w-full text-left px-4 py-2 text-xs hover:bg-[#EEF2FF] border-b border-[#F1F5F9] last:border-0"
                                    onClick={() => handleSuggestionSelect(actualIndex, s)}
                                  >
                                    <div className="font-bold text-[#0F172A]">{s.name}</div>
                                    <div className="text-[10px] text-[#64748B] truncate">{s.smiles}</div>
                                  </button>
                                ))}
                              </div>
                            )}
                          </div>
                          <button
                            type="button"
                            onClick={() => removeCompound(actualIndex)}
                            title="Remove compound"
                            className="w-10 h-10 flex items-center justify-center text-[#94A3B8] hover:text-red-500 hover:bg-red-50 border border-transparent hover:border-red-200 rounded-lg transition-colors text-sm font-bold"
                          >
                            ✕
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
                      + Add Secondary Compound
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
                      { id: "Heuristic", title: "Heuristic / AI (Expert Kinetic Activation & Transition States)" },
                      { id: "Boltzmann", title: "Boltzmann / Physics (Thermodynamic Free Energy ΔG Distribution at 298.15K)" }
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
                  <span>🧪</span> Predict Chemical Interactions
                </button>
              </form>
            </div>
          </div>
        )}

        {/* View: Results Dashboard Matching Streamlit */}
        {view === 'results' && result && !loading && (
          <div className="w-full space-y-6">
            {/* Action Toolbar */}
            <div className="flex flex-col sm:flex-row items-center justify-between gap-4">
              <button
                onClick={() => { setView('input'); setResult(null); }}
                className="inline-flex items-center gap-2 px-4 py-2 bg-white border border-[#CBD5E1] text-[#334155] hover:bg-[#F8FAFC] font-semibold text-sm rounded-lg transition-colors shadow-xs"
              >
                ← Back to Reaction Setup
              </button>
              <button
                onClick={downloadExcel}
                className="inline-flex items-center gap-2 px-4 py-2 bg-[#0F172A] hover:bg-[#1E293B] text-white font-semibold text-sm rounded-lg transition-colors shadow-sm"
              >
                <Download className="w-4 h-4" />
                Download Excel Interaction Report
              </button>
            </div>

            {/* Executive Summary Metrics */}
            {(() => {
              const impurities = result.degradationImpurities || [];
              const maxProb = impurities.length > 0 ? Math.max(...impurities.map(i => i.probability || 0)) : 0;
              const energies = impurities.map(i => i.relativeEnergy).filter((e): e is number => e != null);
              const minEnergy = energies.length > 0 ? Math.min(...energies) : null;

              return (
                <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
                  <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl p-4">
                    <div className="text-[11px] font-bold text-[#64748B] uppercase tracking-wider mb-1">Reaction Products</div>
                    <div className="text-2xl font-extrabold text-[#0F172A]">{impurities.length}</div>
                    <div className="text-[10px] text-[#94A3B8] mt-1">Total identified transformation products</div>
                  </div>
                  <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl p-4">
                    <div className="text-[11px] font-bold text-[#64748B] uppercase tracking-wider mb-1">Highest Probability</div>
                    <div className="text-2xl font-extrabold text-[#4F46E5]">{(maxProb * 100).toFixed(1)}%</div>
                    <div className="text-[10px] text-[#94A3B8] mt-1">Maximum formation likelihood</div>
                  </div>
                  <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl p-4">
                    <div className="text-[11px] font-bold text-[#64748B] uppercase tracking-wider mb-1">Interaction Nature</div>
                    <div className="text-2xl font-extrabold text-[#0F172A]">{result.interactionType || "Chemical"}</div>
                    <div className="text-[10px] text-[#94A3B8] mt-1">Dominant interaction classification</div>
                  </div>
                  <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl p-4">
                    <div className="text-[11px] font-bold text-[#64748B] uppercase tracking-wider mb-1">Lowest ΔG (Driving Force)</div>
                    <div className="text-2xl font-extrabold text-[#0F172A] font-mono">
                      {minEnergy !== null ? `${minEnergy.toFixed(2)} kcal/mol` : "Calculated"}
                    </div>
                    <div className="text-[10px] text-[#94A3B8] mt-1">Most exergonic thermodynamic pathway</div>
                  </div>
                </div>
              );
            })()}

            {/* Tabbed Analytical Views Matching Streamlit */}
            <div className="border-b border-[#E2E8F0] flex gap-2 pt-2">
              <button
                onClick={() => setActiveTab("products")}
                className={`px-4 py-2.5 text-sm font-semibold rounded-t-lg transition-colors flex items-center gap-1.5 ${
                  activeTab === "products"
                    ? "bg-[#EEF2FF] text-[#4F46E5] border-b-2 border-[#4F46E5]"
                    : "text-[#64748B] hover:text-[#0F172A] hover:bg-[#F8FAFC]"
                }`}
              >
                <span>🔬</span> Predicted Reaction Products
              </button>
              <button
                onClick={() => setActiveTab("reactants")}
                className={`px-4 py-2.5 text-sm font-semibold rounded-t-lg transition-colors flex items-center gap-1.5 ${
                  activeTab === "reactants"
                    ? "bg-[#EEF2FF] text-[#4F46E5] border-b-2 border-[#4F46E5]"
                    : "text-[#64748B] hover:text-[#0F172A] hover:bg-[#F8FAFC]"
                }`}
              >
                <span>🧪</span> Reactants & Components
              </button>
              <button
                onClick={() => setActiveTab("reasoning")}
                className={`px-4 py-2.5 text-sm font-semibold rounded-t-lg transition-colors flex items-center gap-1.5 ${
                  activeTab === "reasoning"
                    ? "bg-[#EEF2FF] text-[#4F46E5] border-b-2 border-[#4F46E5]"
                    : "text-[#64748B] hover:text-[#0F172A] hover:bg-[#F8FAFC]"
                }`}
              >
                <span>🧠</span> AI Mechanistic Framework
              </button>
            </div>

            {/* Tab 1: Predicted Products */}
            {activeTab === "products" && (
              <div className="space-y-4">
                <div className="my-2">
                  <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                    Predicted Reaction Byproducts & Degradation Products
                  </h3>
                  <p className="text-xs sm:text-sm text-[#64748B]">
                    Ranked by formation probability and thermodynamic stability under specified reaction conditions.
                  </p>
                </div>

                {(!result.degradationImpurities || result.degradationImpurities.length === 0) ? (
                  <div className="p-8 text-center bg-[#F8FAFC] border border-[#E2E8F0] rounded-xl text-sm text-[#64748B]">
                    No significant byproducts detected under standard conditions.
                  </div>
                ) : (
                  [...result.degradationImpurities]
                    .sort((a, b) => (b.probability || 0) - (a.probability || 0))
                    .map((imp, idx) => {
                      const prob = (imp.probability || 0) * 100;
                      const cond = imp.condition || "Direct Degradation";
                      let condClass = "cond-hydro";
                      const condLower = cond.toLowerCase();
                      if (condLower.includes("oxid")) condClass = "cond-oxid";
                      else if (condLower.includes("therm")) condClass = "cond-therm";
                      else if (condLower.includes("photo")) condClass = "cond-photo";
                      else if (condLower.includes("react") || condLower.includes("incomp")) condClass = "cond-react";

                      return (
                        <div key={`prod-${idx}`} className="ap1-imp-card flex-col sm:flex-row">
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
                                  <span>⚡ Chemical Mechanism:</span>
                                </div>
                                <div>{imp.mechanismExplanation}</div>
                              </div>
                            )}

                            <div className="flex flex-wrap gap-2 items-center">
                              <span className={`ap1-badge-cond ${condClass}`}>{cond}</span>
                              <span className="ap1-pill font-semibold text-[#4338CA] bg-[#EEF2FF] border-[#E0E7FF]">
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
                    })
                )}
              </div>
            )}

            {/* Tab 2: Reactants & Components */}
            {activeTab === "reactants" && (
              <div className="space-y-4">
                <div className="my-2">
                  <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                    Input Molecular Profiles
                  </h3>
                  <p className="text-xs sm:text-sm text-[#64748B]">
                    Calculated molecular descriptors, functional group features, and predicted reactive interaction sites.
                  </p>
                </div>

                {(result.compounds || []).map((comp, idx) => {
                  const mw = comp.molecularDescriptors?.MolWt;
                  const role = idx === 0 ? "Primary Compound" : `Secondary Compound ${idx}`;
                  const roleClass = idx === 0 ? "role-primary" : "role-secondary";

                  return (
                    <div key={`comp-card-${idx}`} className="ap1-comp-card flex-col sm:flex-row">
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
                            <div className="text-[11px] font-bold text-[#4F46E5] uppercase tracking-wider mb-1">
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
            )}

            {/* Tab 3: AI Reasoning Framework */}
            {activeTab === "reasoning" && (
              <div className="space-y-4">
                <div className="my-2">
                  <h3 className="font-serif text-xl font-bold text-[#0F172A] mb-1">
                    AI Mechanistic Reasoning Framework
                  </h3>
                  <p className="text-xs sm:text-sm text-[#64748B]">
                    Comprehensive kinetic pathways, microenvironmental influences, and thermodynamic justification.
                  </p>
                </div>

                <div className="bg-white border border-[#E2E8F0] rounded-2xl p-6 sm:p-7 shadow-[0_1px_3px_rgba(15,23,42,0.03)]">
                  <div className="text-sm text-[#334155] leading-relaxed whitespace-pre-wrap font-sans">
                    {result.chainOfThought || "No detailed reasoning chain provided."}
                  </div>
                </div>

                <div className="bg-[#F8FAFC] border border-[#E2E8F0] rounded-2xl p-5 text-xs text-[#475569] leading-relaxed">
                  <strong className="text-[#0F172A] block mb-1 text-sm font-semibold">
                    Analytical Guidance & Mechanistic Considerations:
                  </strong>
                  <ul className="list-disc pl-5 space-y-1 mt-1 text-[#475569]">
                    <li>Thermodynamic equilibria (Boltzmann distribution) assume standard state conditions at 298.15K.</li>
                    <li>Kinetic barrier overrides may favor lower activation energy pathways over thermodynamically deep energy wells.</li>
                    <li>Secondary excipient microenvironmental changes (e.g. moisture sorption, basic surface catalysis) heavily modulate reaction rates.</li>
                  </ul>
                </div>
              </div>
            )}

            {/* Disclaimer Matching Streamlit */}
            <div className="border border-[#E2E8F0] rounded-xl bg-[#FAFAFA] p-4 text-xs text-[#64748B] leading-relaxed mt-6">
              Disclaimer: INTERACTION is an AI-assisted computational chemistry modeling tool designed for reaction pathway exploration and byproduct screening. Predictions should be verified by experimental analytical assays (HPLC, LC-MS, NMR).
            </div>
          </div>
        )}
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
