import { GoogleGenAI, Type } from "@google/genai";
import { parse as parsePartial } from 'partial-json';
import { validateSmiles, validateSmarts, MolecularDescriptors } from "./rdkit";

const getApiKey = () => {
  // Try both standard and VITE_ prefixed variables
  const key = process.env.GEMINI_API_KEY || process.env.VITE_GEMINI_API_KEY;
  
  if (!key || key === "INVALID_OR_MISSING_KEY") {
    if (typeof window !== "undefined") {
      console.error("❌ GEMINI_API_KEY is missing.");
      console.log("1. AI Studio: Click the gear icon (Settings) -> Secrets -> Add GEMINI_API_KEY.");
      console.log("2. Vercel: Add VITE_GEMINI_API_KEY to your Project Settings and redeploy.");
    }
    return "MISSING_KEY";
  }
  return key;
};

const ai = new GoogleGenAI({ apiKey: getApiKey() });

export interface CompoundInfo {
  name: string;
  smiles: string;
  features: string[];
  interactionSites: string[];
  molecularDescriptors?: MolecularDescriptors;
}

export interface PredictionResult {
  chainOfThought: string;
  compounds: CompoundInfo[];
  interactionType: "Physical" | "Chemical" | "None";
  mechanism: string;
  degradationImpurities: {
    iupacName: string;
    smiles: string;
    structureDescription: string;
    origin: string;
    probability: number; // Primary probability for ranking
    probabilityHeuristic?: number;
    probabilityBoltzmann?: number;
    relativeEnergy?: number;
    condition: "Oxidation" | "Acidic Hydrolysis" | "Basic Hydrolysis" | "Photodegradation" | "Thermal Degradation";
    source: "Stress degradation" | "Interaction with other compound";
    mechanismExplanation: string;
    molecularDescriptors?: MolecularDescriptors;
  }[];
}

export type PredictionMethod = "Boltzmann" | "Heuristic" | "Both";

export type InputType = "Name" | "SMILES" | "SMARTS" | "InChI";

export class AnalysisError extends Error {
  constructor(public message: string, public type: string) {
    super(message);
    this.name = "AnalysisError";
  }
}

export interface CompoundInput {
  value: string;
  type: InputType;
  descriptors?: MolecularDescriptors | null;
  originalName?: string;
}

export async function predictInteraction(
  inputs: CompoundInput[],
  method: PredictionMethod = "Heuristic",
  onChunk?: (partialResult: Partial<PredictionResult>) => void
): Promise<PredictionResult> {
  // Pre-validation with RDKit
  for (const input of inputs) {
    if (input.type === "SMILES") {
      const smiles = input.value.trim();
      const validation = await validateSmiles(smiles);
      if (!validation.isValid) {
        throw new AnalysisError(validation.error || "Invalid chemical structure", "INVALID_SMILES");
      }
      // Use canonical SMILES for better model performance
      if (validation.canonicalSmiles) {
        input.value = validation.canonicalSmiles;
      }
    } else if (input.type === "SMARTS") {
      const smarts = input.value.trim();
      const validation = await validateSmarts(smarts);
      if (!validation.isValid) {
        throw new AnalysisError(validation.error || "Invalid SMARTS pattern", "INVALID_SMARTS");
      }
    }
  }

  const compoundsInfo = inputs
    .map((input, i) => {
      let desc = "";
      if (input.descriptors) {
         desc = ` [Calculated Specs: MolWt: ${input.descriptors.MolWt?.toFixed(2) || 'N/A'}, LogP: ${input.descriptors.MolLogP?.toFixed(2) || 'N/A'}, TPSA: ${input.descriptors.TPSA?.toFixed(2) || 'N/A'}, Rotatable Bonds: ${input.descriptors.NumRotatableBonds ?? 'N/A'}]`;
      }
      const compoundTargetName = input.originalName || (input.type === "Name" ? input.value : `Compound ${i + 1}`);
      const rawType = input.type === "SMILES" ? `SMILES Structure Data: ${input.value}` : `"${input.value}" (provided as Name)`;
      const constraint = ` The user specifically named this compound "${compoundTargetName}". You MUST strictly fill out the 'name' field using exactly "${compoundTargetName}"... DO NOT under any circumstances output 'Compound ${i + 1}' or its IUPAC name for the Input Compound name.`;
      return `Compound ${i + 1}: ${rawType}.${constraint}${desc}`;
    })
    .join("\n    ");

  let probabilityInstruction = "";
  if (method === "Boltzmann") {
    probabilityInstruction = "Probability of formation based on Boltzmann distribution at 298.15K. You MUST also provide the estimated relative formation energy (relativeEnergy) in kcal/mol.";
  } else if (method === "Heuristic") {
    probabilityInstruction = "Probability of formation based on chemical stability principles and heuristic reasoning.";
  } else if (method === "Both") {
    probabilityInstruction = "Provide BOTH 'probabilityHeuristic' (based on expert reasoning) and 'probabilityBoltzmann' (based on thermodynamic ΔG at 298.15K). You MUST also provide 'relativeEnergy' in kcal/mol. The main 'probability' field should match 'probabilityBoltzmann' for ranking purposes.";
  }

  const impurityProperties: any = {
    iupacName: { type: Type.STRING, description: "The IUPAC name or common name of the NEW DEGRADATION PRODUCT. DO NOT just output the parent compound name. You must name the new impurity generated." },
    smiles: { 
      type: Type.STRING, 
      description: "SMILES string of the newly formed degradant. CRITICAL WARNING: You must mathematically ensure standard valence rules are obeyed. Do not attach 5 bonds to Carbon. RDKit will fail to parse this if valences are exceeded."
    },
    structureDescription: { type: Type.STRING },
    origin: { type: Type.STRING, description: "Which specific compound(s) this impurity originated from. E.g. 'Compound 1 and Compound 2'" },
    probability: { 
      type: Type.NUMBER, 
      description: "Primary probability of formation as a decimal between 0.0 and 1.0 (e.g., 0.85 for 85%). Used for ranking." 
    },
    condition: { 
      type: Type.STRING, 
      enum: ["Oxidation", "Acidic Hydrolysis", "Basic Hydrolysis", "Photodegradation", "Thermal Degradation"] 
    },
    source: { 
      type: Type.STRING, 
      enum: ["Stress degradation", "Interaction with other compound"] 
    },
    mechanismExplanation: {
      type: Type.STRING,
      description: "Brief explanation of the interaction mechanism (e.g., pH change, oxidation, complexation, adsorption, precipitation). Mention effects on stability or release kinetics."
    }
  };

  const requiredImpurityFields = ["iupacName", "smiles", "structureDescription", "origin", "probability", "condition", "source", "mechanismExplanation"];

  if (method === "Boltzmann" || method === "Both") {
    impurityProperties.relativeEnergy = { type: Type.NUMBER, description: "Relative formation energy in kcal/mol" };
    requiredImpurityFields.push("relativeEnergy");
  }

  if (method === "Both") {
    impurityProperties.probabilityHeuristic = { 
      type: Type.NUMBER, 
      description: "Heuristic-based probability as a decimal between 0.0 and 1.0" 
    };
    impurityProperties.probabilityBoltzmann = { 
      type: Type.NUMBER, 
      description: "Boltzmann-based probability as a decimal between 0.0 and 1.0" 
    };
    requiredImpurityFields.push("probabilityHeuristic", "probabilityBoltzmann");
  }

  try {
    const responseStream = await ai.models.generateContentStream({
      model: "gemini-2.5-pro",
      contents: `Predict and evaluate the degradation of Compound 1 in the following mixture using the ${method === "Both" ? "Heuristic AND Boltzmann" : method}-based approach:\n${compoundsInfo}`,
      config: {
        temperature: 0.1,
        systemInstruction: `You are a professional pharmaceutical degradation evaluator. 
        Your task is to predict the degradation of Compound 1 using the following analytical framework${method === "Both" ? "s" : ""}:
        ${method === "Heuristic" || method === "Both" ? "\n        1. HEURISTIC ANALYSIS: Based on expert chemical reasoning, reactive site identification, and known reaction kinetics." : ""}
        ${method === "Boltzmann" || method === "Both" ? `\n        ${method === "Both" ? "2." : "1."} BOLTZMANN ANALYSIS: Based on thermodynamic stability and calculated relative formation energy (ΔG) at 298.15K.` : ""}
        
        ${method === "Both" ? "When \"Both\" is selected, you must perform these two analyses independently for each predicted impurity to provide a comparative perspective." : ""}
        
        Evaluate the degradation of Compound 1 due to:
        1. Stress degradation of Compound 1.
        2. Interactions between Compound 1 and any other provided compounds (Compounds 2-5).
        
        You MUST evaluate degradation under these specific conditions:
        - Oxidation
        - Acidic Hydrolysis
        - Basic Hydrolysis
        - Photodegradation
        - Thermal Degradation
        
        First, identify the chemical structures correctly for ALL provided compounds.
        For each compound, provide:
        - Identified name (If the user explicitly provided a name, you MUST echo their exact original name back to them. DO NOT rename it to IUPAC or another common name).
        - SMILES string (MUST be a valid, standard, canonical SMILES string compatible with RDKit and PubChem).
        - List of key structural features.
        - List of specific "Interaction Sites" likely to be involved in degradation.
        
        Predict exactly 10 possible degradation impurities or interaction products derived from Compound 1.
        For each product, you MUST specify:
        - Whether it forms from "Stress degradation" or "Interaction with other compound".
        - Which specific condition it forms under.
        - IUPAC name of the NEW DEGRADANT (do NOT just repeat the starting material's name. You must identify the unique name of the resulting product).
        - SMILES string of the new degradant. CRITICAL: This MUST be a valid, canonical SMILES string. You MUST implicitly verify that standard valences are not exceeded (e.g. Carbon max 4 bonds, Oxygen max 2 bonds, Nitrogen max 3 or 4 if charged) so that RDKit can successfully parse and render it. Invalid SMILES will break the UI renderer. If uncertain, simplify the resulting structure to ensure valence validity.
        - A brief explanation of the underlying mechanism (mechanismExplanation), including specific phenomena like pH changes, oxidation, complexation, adsorption, precipitation, or effects on stability and release kinetics.
        - ${probabilityInstruction}
        
        IMPORTANT: Probabilities MUST be realistic estimates between 0.01 and 0.99. DO NOT return 0.0 unless the impurity is chemically impossible.
        
        Rank the products by their calculated Boltzmann probability (if available) or general probability.`,
        responseMimeType: "application/json",
        responseSchema: {
          type: Type.OBJECT,
          properties: {
            chainOfThought: {
              type: Type.STRING,
              description: `Perform your step-by-step chemical reasoning, mechanism formulation${method !== "Heuristic" ? ", and energy estimation" : ""} here BEFORE outputting the final compounds.`
            },
            compounds: {
              type: Type.ARRAY,
              items: {
                type: Type.OBJECT,
                properties: {
                  name: { 
                    type: Type.STRING,
                    description: "You MUST strictly echo the user's original name exactly as it was provided. Do not convert the starting materials into IUPAC names."
                  },
                  smiles: { type: Type.STRING },
                  features: { type: Type.ARRAY, items: { type: Type.STRING } },
                  interactionSites: { 
                    type: Type.ARRAY, 
                    items: { type: Type.STRING },
                    description: "Specific chemical sites/groups likely to interact or degrade"
                  }
                },
                required: ["name", "smiles", "features", "interactionSites"]
              }
            },
            interactionType: { type: Type.STRING, enum: ["Physical", "Chemical", "None"] },
            mechanism: { type: Type.STRING },
            degradationImpurities: {
              type: Type.ARRAY,
              items: {
                type: Type.OBJECT,
                properties: impurityProperties,
                required: requiredImpurityFields
              }
            }
          },
          required: ["chainOfThought", "compounds", "interactionType", "mechanism", "degradationImpurities"]
        }
      }
    });

    let fullText = "";
    
    for await (const chunk of responseStream) {
      if (chunk.text) {
        fullText += chunk.text;
        if (onChunk) {
          try {
            // Use partial-json to parse the incomplete JSON string
            const partial = parsePartial(fullText);
            if (partial) {
              onChunk(partial as Partial<PredictionResult>);
            }
          } catch (e) {
            // Ignore partial parse errors for intermediate chunks
          }
        }
      }
    }

    if (!fullText) {
      throw new AnalysisError("The model failed to generate a response. Please try again.", "EMPTY_RESPONSE");
    }

    try {
      return JSON.parse(fullText);
    } catch (e) {
      console.error("JSON Parse Error:", fullText);
      throw new AnalysisError("The model generated an invalid chemical report. This can happen with very complex structures. Please try again.", "INVALID_JSON");
    }
  } catch (error: any) {
    if (error instanceof AnalysisError) throw error;

    console.error("Gemini API Error:", error);
    
    // Handle specific Gemini error cases
    if (error.message?.includes("SAFETY")) {
      throw new AnalysisError("The input contains content that triggered safety filters. Please check your compound names or structures.", "SAFETY_TRIGGERED");
    }

    // Check for rate limit (429) or quota errors
    const isRateLimit = error.message?.includes("quota") || 
                        error.message?.includes("429") || 
                        error.status === 429 ||
                        error.message?.toLowerCase().includes("rate limit");

    if (isRateLimit) {
      throw new AnalysisError("The daily request quota for the Gemini API has been reached or you are sending requests too quickly. Please wait 60 seconds and try again.", "QUOTA_EXCEEDED");
    }

    if (error.message?.includes("network") || error.message?.includes("fetch")) {
      throw new AnalysisError("A network connection issue was detected. Please check your internet connection and verify that the API is reachable.", "CONNECTION_ERROR");
    }
    
    if (error.message?.includes("overloaded") || error.status === 503) {
      throw new AnalysisError("The AI engine is currently experiencing high volume and is temporarily overloaded. Please try your request again in a few moments.", "MODEL_OVERLOADED");
    }

    if (error.message?.includes("API key not valid")) {
      throw new AnalysisError("The configured Gemini API key is invalid or has been revoked. Please verify your credentials in the Settings menu.", "CONFIG_ERROR");
    }

    // Default error with more context
    const errorMessage = error.message || "An analytical failure occurred while processing the molecular structures.";
    throw new AnalysisError(`${errorMessage}`, "UNKNOWN_ERROR");
  }
}
