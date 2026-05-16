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
  const [selectedMethods, setSelectedMethods] = useState<Set<"Heuristic" | "Boltzmann">>(new Set(["Heuristic"]));
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

  const toggleMethod = (method: "Heuristic" | "Boltzmann") => {
    const newMethods = new Set(selectedMethods);
    if (newMethods.has(method)) {
      newMethods.delete(method);
    } else {
      newMethods.add(method);
    }
    setSelectedMethods(newMethods);
  };

  const predictionMethod: PredictionMethod = 
    selectedMethods.size === 2 ? "Both" : 
    selectedMethods.has("Boltzmann") ? "Boltzmann" : 
    selectedMethods.has("Heuristic") ? "Heuristic" : 
    "Heuristic";

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
      rows.push(["A-Pi1 REPORT"]);
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
        rows.push(["PREDICTED DEGRADATION IMPURITIES"]);
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

      XLSX.utils.book_append_sheet(wb, ws, "Stability Report");
      XLSX.writeFile(wb, `A-Pi1_Report_${new Date().toISOString().split('T')[0]}.xlsx`);
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
      <div className="min-h-screen bg-background flex flex-col">
      {/* Header */}
      <header className="border-b bg-white sticky top-0 z-10">
        <div className="max-w-7xl mx-auto px-4 h-16 flex items-center justify-between">
          <div className="flex items-center gap-2">
            <h1 className="text-xl font-bold tracking-tight font-serif">A-Pi<span className="text-2xl">1</span></h1>
          </div>
          <div className="flex items-center gap-4">
            <div className="hidden md:flex items-center gap-3 px-3 py-1 bg-slate-50 rounded-full border border-slate-100">
              <div className="flex items-center gap-1.5">
                <Database className="w-3 h-3 text-indigo-500" />
                <span className="text-[10px] font-medium text-slate-500">{dbStats.compounds} Compounds</span>
              </div>
              <Separator orientation="vertical" className="h-3 bg-slate-200" />
              <div className="flex items-center gap-1.5">
                <History className="w-3 h-3 text-emerald-500" />
                <span className="text-[10px] font-medium text-slate-500">{dbStats.predictions} Predictions</span>
              </div>
          </div>
        </div>
      </div>
    </header>


      <main className="flex-1 max-w-7xl mx-auto px-4 py-8 w-full relative">
        {view === 'loading' && (
          <div className="fixed inset-0 z-50 flex flex-col items-center justify-center bg-white">
            <div className="relative w-32 h-32 flex items-center justify-center mb-8">
              {/* Nucleus */}
              <div className="w-8 h-8 bg-indigo-600 rounded-full shadow-[0_0_30px_#4f46e5] z-10" />
              
              {/* Orbit Path */}
              <div className="absolute w-24 h-24 border-2 border-slate-100 rounded-full" />
              
              {/* Electron Container (Rotates) */}
              <motion.div 
                className="absolute w-24 h-24"
                animate={{ rotate: 360 }}
                transition={{ duration: 1.5, repeat: Infinity, ease: "linear" }}
              >
                {/* The Electron */}
                <div className="absolute top-0 left-1/2 -translate-x-1/2 -translate-y-1/2 w-4 h-4 bg-indigo-400 rounded-full shadow-[0_0_15px_#818cf8]" />
              </motion.div>
            </div>
            
            <div className="text-center space-y-2">
              <h2 className="text-2xl font-bold text-slate-800">Generating Possible Degradation Products...</h2>
            </div>
          </div>
        )}
        
        {view === 'input' && (
          <div className="w-full max-w-3xl mx-auto space-y-6">
          <Card className="shadow-sm border-slate-200">
            <CardHeader>
              <CardTitle className="text-lg font-serif">Input</CardTitle>
              <CardDescription>Predict interactions and degradation pathways between compounds or analyze intrinsic stability.</CardDescription>
              {error && (
                <div className="pt-4">
                  <Alert variant="destructive" className="border-red-200 bg-red-50">
                    <div className="flex gap-3">
                      <div className="mt-0.5">
                        {(error.type === "QUOTA_EXCEEDED" || error.type === "MODEL_OVERLOADED") && <Clock className="h-5 w-5 text-red-600 animate-pulse" />}
                        {error.type === "SAFETY_TRIGGERED" && <ShieldAlert className="h-5 w-5 text-red-600" />}
                        {error.type === "CONNECTION_ERROR" && <WifiOff className="h-5 w-5 text-red-600" />}
                        {(error.type === "INVALID_SMILES" || error.type === "INVALID_SMARTS" || error.type === "INVALID_JSON") && <AlertCircle className="h-5 w-5 text-red-600" />}
                        {error.type === "CONFIG_ERROR" && <Lock className="h-5 w-5 text-red-600" />}
                        {error.type === "UNKNOWN_ERROR" && <AlertTriangle className="h-5 w-5 text-red-600" />}
                      </div>
                      <div className="space-y-1">
                        <AlertTitle className="text-red-800 font-bold">
                          {error.type === "QUOTA_EXCEEDED" ? "API Quota Exceeded" : 
                           error.type === "MODEL_OVERLOADED" ? "Model Temporarily Overloaded" :
                           error.type === "SAFETY_TRIGGERED" ? "Safety Filter Logic Engaged" :
                           error.type === "INVALID_SMILES" || error.type === "INVALID_SMARTS" ? "Structural Format Error" :
                           error.type === "CONNECTION_ERROR" ? "Network Communication Failure" :
                           error.type === "CONFIG_ERROR" ? "Configuration Credential Error" :
                           error.type === "INVALID_JSON" ? "Structure Interpretation Failure" :
                           "Analytical Processing Error"}
                        </AlertTitle>
                        <AlertDescription className="text-red-700">
                          {error.message}
                        </AlertDescription>
                        {(error.type === "QUOTA_EXCEEDED" || error.type === "MODEL_OVERLOADED" || error.type === "CONNECTION_ERROR") && (
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
            </CardHeader>
            <CardContent>
              <form onSubmit={handlePredict} className="space-y-6" autoComplete="off">
                <div className="space-y-6">
                  {/* Primary Compound Tile */}
                  <div className="space-y-3 p-4 bg-[#f5f7ff] rounded-xl border border-indigo-100 relative group">
                    <div className="flex items-center justify-between">
                      <Label htmlFor={`compound-0`} className="font-bold text-indigo-900 flex items-center gap-2">
                        <div className="w-2 h-2 rounded-full bg-indigo-500" />
                        Primary Compound
                      </Label>
                    </div>
                    <Tabs value={compounds[0].type} onValueChange={(v) => updateCompound(0, "type", v as InputType)} className="w-full">
                      <TabsList className="grid grid-cols-4 w-full h-7">
                        <TabsTrigger value="Name" className="text-[9px]">Name</TabsTrigger>
                        <TabsTrigger value="SMILES" className="text-[9px]">SMILES</TabsTrigger>
                        <TabsTrigger value="SMARTS" className="text-[9px]">SMARTS</TabsTrigger>
                        <TabsTrigger value="InChI" className="text-[9px]">InChI</TabsTrigger>
                      </TabsList>
                    </Tabs>
                    <div className="relative">
                      <Input 
                        id={`compound-0`} 
                        placeholder={compounds[0].type === "Name" ? "e.g., Aspirin" : `Enter ${compounds[0].type}...`}
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
                        className="bg-white h-9 text-sm border-indigo-100 focus-visible:ring-indigo-500"
                      />
                      {showSuggestions === 0 && suggestions.length > 0 && (
                        <div className="absolute z-50 w-full mt-1 bg-white border border-slate-200 rounded-lg shadow-xl overflow-hidden max-h-48 overflow-y-auto">
                          {suggestions.map((s, i) => (
                            <button
                              key={i}
                              type="button"
                              className="w-full text-left px-4 py-2 text-xs hover:bg-indigo-50 border-b border-slate-50 last:border-0"
                              onClick={() => handleSuggestionSelect(0, s)}
                            >
                              <div className="font-bold text-slate-700">{s.name}</div>
                              <div className="text-[10px] text-slate-400 truncate">{s.smiles}</div>
                            </button>
                          ))}
                        </div>
                      )}
                    </div>
                  </div>

                  {/* Secondary Compounds Tile */}
                  <div className="space-y-4 p-4 bg-[#f8fafc] rounded-xl border border-slate-100">
                    <div className="flex items-center justify-between">
                      <Label className="font-semibold text-slate-700">Secondary Compounds</Label>
                      <Badge variant="outline" className="text-[9px] font-mono text-slate-400">{compounds.length - 1} Added</Badge>
                    </div>
                    
                    <div className="space-y-4">
                      {compounds.slice(1).map((c, idx) => {
                        const actualIndex = idx + 1;
                        return (
                          <div key={`input-compound-${actualIndex}`} className="space-y-3 relative group">
                            <div className="flex items-center justify-between">
                              <span className="text-[10px] font-bold text-slate-400 uppercase tracking-wider">Secondary Compound {idx + 1}</span>
                              <Button 
                                type="button" 
                                variant="ghost" 
                                size="icon" 
                                className="h-5 w-5 text-slate-300 hover:text-red-500"
                                onClick={() => removeCompound(actualIndex)}
                              >
                                <Trash2 className="h-3 w-3" />
                              </Button>
                            </div>
                            <Tabs value={c.type} onValueChange={(v) => updateCompound(actualIndex, "type", v as InputType)} className="w-full">
                              <TabsList className="grid grid-cols-4 w-full h-6">
                                <TabsTrigger value="Name" className="text-[8px]">Name</TabsTrigger>
                                <TabsTrigger value="SMILES" className="text-[8px]">SMILES</TabsTrigger>
                                <TabsTrigger value="SMARTS" className="text-[8px]">SMARTS</TabsTrigger>
                                <TabsTrigger value="InChI" className="text-[8px]">InChI</TabsTrigger>
                              </TabsList>
                            </Tabs>
                            <div className="relative">
                              <Input 
                                placeholder={c.type === "Name" ? "e.g., Lactose" : `Enter ${c.type}...`}
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
                                className="bg-white h-8 text-xs"
                              />
                              {showSuggestions === actualIndex && suggestions.length > 0 && (
                                <div className="absolute z-50 w-full mt-1 bg-white border border-slate-200 rounded-lg shadow-xl overflow-hidden max-h-48 overflow-y-auto">
                                  {suggestions.map((s, i) => (
                                    <button
                                      key={i}
                                      type="button"
                                      className="w-full text-left px-4 py-2 text-xs hover:bg-indigo-50 border-b border-slate-50 last:border-0"
                                      onClick={() => handleSuggestionSelect(actualIndex, s)}
                                    >
                                      <div className="font-bold text-slate-700">{s.name}</div>
                                      <div className="text-[10px] text-slate-400 truncate">{s.smiles}</div>
                                    </button>
                                  ))}
                                </div>
                              )}
                            </div>
                          </div>
                        );
                      })}
                    </div>

                    {compounds.length < 5 && (
                      <Button 
                        type="button" 
                        variant="outline" 
                        className="w-full border-dashed border-slate-300 text-slate-500 hover:text-indigo-600 hover:border-indigo-300 h-9 text-xs"
                        onClick={addCompound}
                      >
                        <Plus className="mr-2 h-3.5 w-3.5" />
                        Add Secondary Compound
                      </Button>
                    )}
                  </div>
                </div>

                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <Label className="text-xs font-bold text-slate-500 uppercase tracking-wider">Prediction Method</Label>
                    <Tooltip>
                      <TooltipTrigger className="text-slate-400 hover:text-indigo-500 transition-colors flex items-center justify-center">
                        <Info className="w-3.5 h-3.5" />
                      </TooltipTrigger>
                        <TooltipContent className="max-w-xs p-3 bg-white border border-slate-200 shadow-xl rounded-xl">
                          <div className="space-y-2">
                            <p className="text-xs font-bold text-indigo-900">Heuristic/AI</p>
                            <p className="text-[10px] text-slate-600 leading-relaxed">Uses expert pharmaceutical reasoning and reaction kinetics to identify reactive sites and rank outcomes.</p>
                            <div className="h-px bg-slate-100" />
                            <p className="text-xs font-bold text-indigo-900">Boltzmann/Physics</p>
                            <p className="text-[10px] text-slate-600 leading-relaxed">Uses thermodynamic stability principles and provides relative formation energy (ΔG) for each predicted impurity.</p>
                          </div>
                        </TooltipContent>
                      </Tooltip>
                  </div>
                  <div className="flex flex-col sm:flex-row gap-2">
                    <button
                      type="button"
                      onClick={() => toggleMethod("Heuristic")}
                      className={`flex-1 flex items-center gap-2.5 p-2.5 rounded-xl border transition-all text-left ${
                        selectedMethods.has("Heuristic")
                          ? "bg-indigo-50 border-indigo-200 ring-1 ring-indigo-200"
                          : "bg-white border-slate-200 hover:border-slate-300"
                      }`}
                    >
                      <div className={`shrink-0 w-4 h-4 rounded-sm border flex items-center justify-center transition-colors ${
                        selectedMethods.has("Heuristic") ? "bg-indigo-600 border-indigo-600" : "bg-white border-slate-300"
                      }`}>
                        {selectedMethods.has("Heuristic") && <div className="w-1.5 h-1.5 bg-white rounded-full" />}
                      </div>
                      <div className="min-w-0">
                        <div className={`text-xs font-bold truncate ${selectedMethods.has("Heuristic") ? "text-indigo-900" : "text-slate-700"}`}>Heuristic/AI</div>
                        <div className="text-[9px] text-slate-500 truncate">Expert Reasoning</div>
                      </div>
                    </button>
                    <button
                      type="button"
                      onClick={() => toggleMethod("Boltzmann")}
                      className={`flex-1 flex items-center gap-2.5 p-2.5 rounded-xl border transition-all text-left ${
                        selectedMethods.has("Boltzmann")
                          ? "bg-indigo-50 border-indigo-200 ring-1 ring-indigo-200"
                          : "bg-white border-slate-200 hover:border-slate-300"
                      }`}
                    >
                      <div className={`shrink-0 w-4 h-4 rounded-sm border flex items-center justify-center transition-colors ${
                        selectedMethods.has("Boltzmann") ? "bg-indigo-600 border-indigo-600" : "bg-white border-slate-300"
                      }`}>
                        {selectedMethods.has("Boltzmann") && <div className="w-1.5 h-1.5 bg-white rounded-full" />}
                      </div>
                      <div className="min-w-0">
                        <div className={`text-xs font-bold truncate ${selectedMethods.has("Boltzmann") ? "text-indigo-900" : "text-slate-700"}`}>Boltzmann/Physics</div>
                        <div className="text-[9px] text-slate-500 truncate">Thermodynamic (ΔG)</div>
                      </div>
                    </button>
                  </div>
                </div>

                <Button 
                  type="submit" 
                  className="w-full h-11" 
                  disabled={loading || compounds.every(c => c.value.trim() === "")}
                >
                  {loading ? (
                    <>
                      <RefreshCw className="mr-2 h-4 w-4 animate-spin" />
                      Analyzing...
                    </>
                  ) : (
                    <>
                      <Search className="mr-2 h-4 w-4" />
                      Predict Interaction and Degradation
                    </>
                  )}
                </Button>
              </form>
            </CardContent>
          </Card>
          </div>
        )}

        {/* View: Results */}
        {view === 'results' && result && !loading && (
          <div className="w-full max-w-5xl mx-auto space-y-6">
            <Button variant="ghost" className="mb-4" onClick={() => {
              setView('input');
              // Clear result when returning to input to prevent flashing old results later
              setResult(null);
            }}>
                <ArrowLeft className="mr-2 h-4 w-4" /> Back to Input
            </Button>
                <div className="flex justify-end">
                  <Button 
                    variant="outline" 
                    size="sm" 
                    onClick={downloadExcel}
                    className="bg-white border-slate-200 text-slate-600 hover:text-indigo-600 hover:border-indigo-200"
                  >
                    <Download className="mr-2 h-4 w-4" />
                    Download Excel Report
                  </Button>
                </div>

                <div ref={reportRef} className="space-y-6 bg-white p-1 rounded-xl">
                  <Card className="border-slate-200 overflow-hidden">
                  <CardHeader className="bg-[#f8fafc] relative pb-10">
                    <div className="absolute top-4 left-1/2 -translate-x-1/2">
                      <Badge variant="outline" className="py-1 px-6 text-xs font-bold uppercase tracking-[0.2em] bg-white border-slate-200 text-[#64748b]">
                        Input
                      </Badge>
                    </div>
                    
                    <div className="grid grid-cols-1 gap-4 pt-8">
                      {(!result.compounds || result.compounds.length === 0 ? compounds.filter(c => c.value.trim() !== "") : result.compounds).map((compound, idx) => (
                        <div key={`compound-${idx}`} className="bg-white border border-slate-200 rounded-xl overflow-hidden flex flex-col md:flex-row items-center md:items-stretch min-h-[192px]">
                          <div className="w-full md:w-48 h-48 bg-slate-50 flex items-center justify-center border-b md:border-b-0 md:border-r border-slate-100 p-4 relative">
                            {('smiles' in compound && compound.smiles) ? (
                              <ChemicalStructure 
                                smiles={compound.smiles as string} 
                                width={160} 
                                height={160} 
                              />
                            ) : (
                              <div className="animate-pulse bg-slate-200 h-24 w-24 rounded-full opacity-50" />
                            )}
                            <div className="absolute top-3 left-3 bg-[#f1f5f9] text-slate-500 text-[10px] font-bold px-2 py-1 rounded-sm">
                              C{idx + 1}
                            </div>
                          </div>
                          <div className="p-6 flex-1 flex flex-col justify-center space-y-3">
                            <div className="flex items-center gap-3">
                              <div className="font-serif text-2xl font-bold text-[#0f172a] leading-tight">
                                {('name' in compound && compound?.name) ? compound.name : <div className="h-6 w-32 bg-slate-200 animate-pulse rounded" />}
                              </div>
                              <Badge variant="outline" className={`text-[10px] font-mono border-slate-100 ${idx === 0 ? 'text-[#4f46e5] bg-[#f5f7ff] border-indigo-100' : 'text-[#94a3b8]'}`}>
                                {idx === 0 ? 'Primary Compound' : `Secondary Compound ${idx}`}
                              </Badge>
                            </div>
                            
                            {('smiles' in compound && compound.smiles) ? (
                               <div className="font-mono text-[11px] text-[#94a3b8] break-all leading-relaxed" title={compound.smiles}>{compound.smiles}</div>
                            ) : (
                               <div className="h-4 w-48 bg-slate-100 animate-pulse rounded" />
                            )}
                            
                            <div className="flex flex-wrap gap-2 pt-2">
                              {('molecularDescriptors' in compound && compound?.molecularDescriptors) && (
                                <Badge variant="outline" className="bg-[#f8fafc] text-[#64748b] border-slate-200 text-[10px] py-0 px-2 font-mono">
                                  Molecular Weight: {compound.molecularDescriptors.MolWt != null ? compound.molecularDescriptors.MolWt.toFixed(2) : "N/A"}
                                </Badge>
                              )}
                              {('features' in compound && compound?.features && compound.features.length > 0) ? compound.features.map((feature, i) => (
                                <Badge key={`feature-${idx}-${i}`} variant="outline" className="bg-[#f8fafc] text-[#475569] border-slate-200 text-[10px] py-0 px-2">
                                  {feature}
                                </Badge>
                              )) : (
                                (!('features' in compound)) && (
                                  <>
                                    <div className="h-4 w-16 bg-slate-100 animate-pulse rounded-full" />
                                    <div className="h-4 w-20 bg-slate-100 animate-pulse rounded-full" />
                                  </>
                                )
                              )}
                            </div>
                            {('interactionSites' in compound && compound?.interactionSites && compound.interactionSites.length > 0) && (
                              <div className="pt-2 space-y-1.5">
                                <div className="text-[10px] font-bold text-[#4f46e5] uppercase tracking-wider">Potential Interaction Sites</div>
                                <div className="flex flex-wrap gap-1.5">
                                  {compound.interactionSites.map((site, i) => (
                                    <Badge key={`site-${idx}-${i}`} variant="secondary" className="bg-[#f5f7ff] text-[#4338ca] border-indigo-100 text-[9px] py-0 px-2">
                                      {site}
                                    </Badge>
                                  ))}
                                </div>
                              </div>
                            )}
                          </div>
                        </div>
                      ))}
                    </div>
                  </CardHeader>
                  <CardContent className="pt-6 space-y-8">
                    <div className="bg-slate-50 border border-slate-200 rounded-lg p-5 min-h-[120px]">
                      <h4 className="text-sm font-semibold flex items-center gap-2 text-[#0f172a] mb-2">
                        AI Reasoning Framework
                      </h4>
                      <div className="text-xs text-slate-600 leading-relaxed max-h-64 overflow-y-auto whitespace-pre-wrap">
                        {result.chainOfThought ? result.chainOfThought : (
                           <div className="animate-pulse space-y-2 mt-2">
                             <div className="h-3 bg-slate-200 rounded w-full"></div>
                             <div className="h-3 bg-slate-200 rounded w-5/6"></div>
                             <div className="h-3 bg-slate-200 rounded w-4/6"></div>
                           </div>
                        )}
                      </div>
                    </div>

                    <>
                      <Separator />
                      <section className="space-y-4">
                        <h4 className="text-sm font-semibold flex items-center gap-2 text-[#0f172a]">
                          Predicted Degradation Impurities
                        </h4>
                        <div className="grid grid-cols-1 gap-6">
                          {(!result.degradationImpurities || result.degradationImpurities.length === 0) ? (
                            Array.from({ length: 3 }).map((_, i) => (
                              <div key={`skeleton-impurity-${i}`} className="bg-white border border-slate-200 rounded-xl overflow-hidden flex flex-col md:flex-row min-h-[256px]">
                                <div className="w-full md:w-64 h-64 bg-slate-50 flex items-center justify-center border-b md:border-b-0 md:border-r border-slate-100 p-4">
                                  <div className="animate-pulse bg-slate-200 h-32 w-32 rounded-full opacity-50" />
                                </div>
                                <div className="p-6 flex-1 space-y-4 w-full">
                                  <div className="flex justify-between items-start">
                                    <div className="space-y-2 w-full">
                                      <div className="h-6 w-48 bg-slate-100 animate-pulse rounded" />
                                      <div className="h-4 w-16 bg-slate-50 animate-pulse rounded" />
                                    </div>
                                    <div className="h-6 w-16 bg-slate-100 animate-pulse rounded" />
                                  </div>
                                  <div className="space-y-2 pt-2">
                                    <div className="h-4 w-full bg-slate-50 animate-pulse rounded" />
                                    <div className="h-4 w-5/6 bg-slate-50 animate-pulse rounded" />
                                  </div>
                                  <div className="bg-[#f8fafc] border border-slate-100 rounded-lg p-3 space-y-2 mt-4">
                                    <div className="h-4 w-24 bg-slate-200 animate-pulse rounded mb-1" />
                                    <div className="h-3 w-full bg-slate-100 animate-pulse rounded" />
                                    <div className="h-3 w-4/5 bg-slate-100 animate-pulse rounded" />
                                  </div>
                                </div>
                              </div>
                            ))
                          ) : (
                            [...(result.degradationImpurities || [])]
                              .sort((a, b) => (b?.probability || 0) - (a?.probability || 0))
                              .map((impurity, i) => (
                              <div key={`impurity-${i}-${impurity?.iupacName || 'unk'}-${(impurity?.smiles || '').slice(0, 10)}`} className="bg-white border border-slate-200 rounded-xl overflow-hidden hover:border-indigo-200 transition-colors flex flex-col md:flex-row">
                                <div className="w-full md:w-64 h-64 bg-white flex items-center justify-center border-b md:border-b-0 md:border-r border-slate-100 p-4 relative">
                                  {impurity?.smiles && (
                                    <ChemicalStructure 
                                      smiles={impurity.smiles} 
                                      width={240} 
                                      height={240} 
                                    />
                                  )}
                                  <div className="absolute top-3 left-3 bg-[#f1f5f9] text-[#64748b] text-[10px] font-bold px-2 py-1 rounded-sm">
                                    #{i + 1}
                                  </div>
                                </div>
                                <div className="p-6 flex-1 space-y-4">
                                  <div className="flex justify-between items-start">
                                    <div className="space-y-1">
                                      <div className="font-bold text-xl text-[#0f172a] break-words">
                                        {impurity.iupacName}
                                      </div>
                                      {impurity.molecularDescriptors && (
                                        <div className="text-[11px] font-mono font-medium text-slate-500 bg-slate-50 border border-slate-100 inline-flex px-2 py-0.5 rounded-sm">
                                          MW: {impurity.molecularDescriptors.MolWt != null ? impurity.molecularDescriptors.MolWt.toFixed(2) : "N/A"}
                                        </div>
                                      )}
                                    </div>
                                    <div className="text-right space-y-1">
                                      <div className="text-lg font-bold text-[#4f46e5]">{impurity.probability != null ? (impurity.probability * 100).toFixed(1) : "N/A"}%</div>
                                      {impurity.probabilityHeuristic != null && impurity.probabilityBoltzmann != null && (
                                        <div className="text-[9px] font-medium text-slate-400 uppercase tracking-tighter">
                                          H: {(impurity.probabilityHeuristic * 100).toFixed(1)}% | B: {(impurity.probabilityBoltzmann * 100).toFixed(1)}%
                                        </div>
                                      )}
                                      {impurity.relativeEnergy != null && (
                                        <div className="text-[10px] font-mono text-slate-400">
                                          ΔG: {impurity.relativeEnergy.toFixed(2)} kcal/mol
                                        </div>
                                      )}
                                    </div>
                                  </div>
                                  <div className="text-sm text-[#475569] leading-relaxed">{impurity.structureDescription}</div>
                                  <div className="bg-[#f8fafc] border border-slate-100 rounded-lg p-3 text-xs text-[#475569] leading-relaxed">
                                    <span className="font-bold text-[#0f172a] block mb-1">Mechanism:</span>
                                    {impurity.mechanismExplanation}
                                  </div>

                                  <div className="flex flex-wrap gap-3 pt-2">
                                    {compounds.length > 1 && (
                                      <div className="text-[10px] bg-[#eef2ff] text-[#4338ca] px-3 py-1 rounded-full font-semibold uppercase tracking-wider">
                                        {impurity.origin}
                                      </div>
                                    )}
                                    <div className="text-[10px] bg-[#fffbeb] text-[#b45309] px-3 py-1 rounded-full font-semibold uppercase tracking-wider">
                                      {impurity.condition}
                                    </div>
                                    <div className="text-[10px] bg-[#ecfdf5] text-[#047857] px-3 py-1 rounded-full font-semibold uppercase tracking-wider">
                                      {impurity.source}
                                    </div>
                                    <div className="text-[10px] bg-[#f1f5f9] text-[#475569] px-3 py-1 rounded-full font-mono break-all" title={impurity.smiles}>
                                      {impurity.smiles}
                                    </div>
                                  </div>
                                </div>
                              </div>
                            ))
                          )}
                        </div>
                      </section>
                    </>
                  </CardContent>
                  <CardFooter className="bg-[#fcfcfc] border-t py-4">
                    <p className="text-[10px] text-[#94a3b8] italic">
                      Disclaimer: This prediction is generated by an AI model and should be used for research purposes only. Always verify with experimental data and professional pharmaceutical consultation.
                    </p>
                  </CardFooter>
                </Card>
              </div>
            </div>
            )}
      </main>

      {/* Footer */}
      <footer className="border-t py-6 bg-white">
        <div className="max-w-7xl mx-auto px-4 flex flex-col md:flex-row justify-between items-center gap-4 text-sm text-slate-500">
          <p>© 2026 A-Pi1 Research Lab. All rights reserved.</p>
          <div className="flex gap-6">
            <a href="#" className="hover:text-foreground transition-colors">Documentation</a>
            <a href="#" className="hover:text-foreground transition-colors">Privacy Policy</a>
            <a href="#" className="hover:text-foreground transition-colors">Contact Support</a>
          </div>
        </div>
      </footer>
    </div>
    </TooltipProvider>
  );
}
