import { GoogleGenAI, Type } from "@google/genai";
import { parse as parsePartial } from 'partial-json';
import { validateSmiles, MolecularDescriptors } from "./rdkit";

const getApiKey = () => {
  // Check runtime globals, URL query parameters, and local storage first
  if (typeof window !== "undefined") {
    const win = window as any;
    if (win.__GEMINI_API_KEY__ && win.__GEMINI_API_KEY__ !== "MISSING_KEY") {
      return win.__GEMINI_API_KEY__;
    }
    const params = new URLSearchParams(window.location.search);
    const urlKey = params.get("gemini_key") || params.get("api_key");
    if (urlKey) {
      try {
        localStorage.setItem("GEMINI_API_KEY", urlKey);
      } catch (_) {}
      return urlKey;
    }
    try {
      const storedKey = localStorage.getItem("GEMINI_API_KEY");
      if (storedKey && storedKey !== "MISSING_KEY") {
        return storedKey;
      }
    } catch (_) {}
  }

  // Try both standard and VITE_ prefixed environment variables
  const key = process.env.GEMINI_API_KEY || process.env.VITE_GEMINI_API_KEY;
  
  if (!key || key === "INVALID_OR_MISSING_KEY") {
    return "MISSING_KEY";
  }
  return key;
};

const getAiClient = () => {
  return new GoogleGenAI({ apiKey: getApiKey() });
};

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

export type InputType = "Name" | "SMILES";

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
    iupacName: { type: Type.STRING, description: "The IUPAC name or common name of the NEW REACTION OR DEGRADATION PRODUCT. DO NOT just output the parent compound name. You must name the new product generated." },
    smiles: { 
      type: Type.STRING, 
      description: "SMILES string of the newly formed product. CRITICAL WARNING: You must mathematically ensure standard valence rules are obeyed. Do not attach 5 bonds to Carbon. RDKit will fail to parse this if valences are exceeded."
    },
    structureDescription: { type: Type.STRING },
    origin: { type: Type.STRING, description: "Which specific compound(s) this product originated from. E.g. 'Compound 1 and Compound 2'" },
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

  const modelsToTry = ["gemini-3.1-flash-lite", "gemini-flash-latest", "gemini-3.8-flash"];
  const maxRetries = 2; // Up to 2 retries per attempt
  let attempt = 0;

  for (let mIdx = 0; mIdx < modelsToTry.length; mIdx++) {
    const activeModel = modelsToTry[mIdx];
    attempt = 0;

    while (attempt <= maxRetries) {
      try {
        const ai = getAiClient();
        const responseStream = await ai.models.generateContentStream({
          model: activeModel,
          contents: `Predict and evaluate the chemical interaction and reaction products of Compound 1 in the following mixture using the ${method === "Both" ? "Heuristic AND Boltzmann" : method}-based approach:\n${compoundsInfo}`,
          config: {
            temperature: 0.1,
            systemInstruction: `You are an expert computational chemist and reaction mechanism evaluator. 
        Your task is to predict the chemical interaction and transformation products of Compound 1 using the following analytical framework${method === "Both" ? "s" : ""}:
        ${method === "Heuristic" || method === "Both" ? "\n        1. HEURISTIC ANALYSIS: Based on expert chemical reasoning, reactive site identification, and known reaction kinetics." : ""}
        ${method === "Boltzmann" || method === "Both" ? `\n        ${method === "Both" ? "2." : "1."} BOLTZMANN ANALYSIS: Based on thermodynamic stability and calculated relative formation energy (ΔG) at 298.15K.` : ""}
        
        ${method === "Both" ? "When \"Both\" is selected, you must perform these two analyses independently for each predicted product to provide a comparative perspective." : ""}
        
        Evaluate the chemical reactivity and transformation of Compound 1 due to:
        1. Direct degradation / intrinsic reactivity of Compound 1.
        2. Chemical interactions between Compound 1 and any other provided co-reactants (Compounds 2-5).
        
        You MUST evaluate reactivity under these specific conditions:
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
        - List of specific "Interaction Sites" likely to be involved in reaction or degradation.
        
        Predict exactly 10 possible reaction byproducts, degradation products, or interaction adducts derived from Compound 1.
        For each product, you MUST specify:
        - Whether it forms from "Direct degradation" or "Interaction with other compound".
        - Which specific condition it forms under.
        - IUPAC name of the NEW PRODUCT (do NOT just repeat the starting material's name. You must identify the unique name of the resulting product).
        - SMILES string of the new product. CRITICAL: This MUST be a valid, canonical SMILES string. You MUST implicitly verify that standard valences are not exceeded (e.g. Carbon max 4 bonds, Oxygen max 2 bonds, Nitrogen max 3 or 4 if charged) so that RDKit can successfully parse and render it. Invalid SMILES will break the UI renderer. If uncertain, simplify the resulting structure to ensure valence validity.
        - A brief explanation of the underlying mechanism (mechanismExplanation), including specific phenomena like pH changes, oxidation, complexation, addition, substitution, rearrangement, or precipitation.
        - ${probabilityInstruction}
        
        IMPORTANT: Probabilities MUST be realistic estimates between 0.01 and 0.99. DO NOT return 0.0 unless the product is chemically impossible.
        
        Rank the products by their calculated Boltzmann probability (if available) or general probability.`,
            responseMimeType: "application/json",
            responseSchema: {
              type: Type.OBJECT,
              properties: {
                chainOfThought: {
                  type: Type.STRING,
                  description: `Perform your step-by-step chemical reasoning, reaction pathway derivation${method !== "Heuristic" ? ", and energy estimation" : ""} here BEFORE outputting the final compounds.`
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
        // If error is AnalysisError we generated, maybe retry if it's transient
        const isInvalidJson = error instanceof AnalysisError && error.type === "INVALID_JSON";
        const isEmpty = error instanceof AnalysisError && error.type === "EMPTY_RESPONSE";

        // Check for rate limit (429) or quota errors
        const isRateLimit = error.message?.includes("quota") || 
                            error.message?.includes("429") || 
                            error.status === 429 ||
                            error.message?.toLowerCase().includes("rate limit") ||
                            error.message?.includes("RESOURCE_EXHAUSTED");

        const isOverloaded = error.message?.includes("overloaded") || 
                             error.message?.includes("high demand") ||
                             error.message?.includes("UNAVAILABLE") ||
                             error.status === 503;

        const isPermissionDenied = error.message?.toLowerCase().includes("permission") ||
                                   error.message?.includes("PERMISSION_DENIED") ||
                                   error.status === 403 ||
                                   error.message?.includes("403") ||
                                   error.message?.includes("unregistered callers") ||
                                   error.message?.includes("not have permission") ||
                                   error.message?.includes("caller does not have permission");

        // If permission denied or other fatal error on this model and we have another model, try next model
        if ((isPermissionDenied || isOverloaded) && mIdx < modelsToTry.length - 1) {
          console.warn(`Model ${activeModel} failed with ${error.message}. Trying next candidate model...`);
          break; // move to next model in modelsToTry
        }
        
        // If we haven't reached max retries and it's a retryable error
        if (attempt < maxRetries && (isRateLimit || isOverloaded || isInvalidJson || isEmpty || error.message?.includes("network") || error.message?.includes("fetch"))) {
          attempt++;
          const delayMs = attempt * 2500 + Math.random() * 1500;
          console.warn(`Transient error encountered on ${activeModel} (${error.message}). Retrying in ${Math.round(delayMs/1000)}s (Attempt ${attempt} of ${maxRetries})...`);
          await new Promise(resolve => setTimeout(resolve, delayMs));
          continue; // Retry the loop
        }

        // If we've exhausted retries or it's a fatal error and we have more models, try next model
        if (mIdx < modelsToTry.length - 1 && (isOverloaded || isRateLimit)) {
          console.warn(`Model ${activeModel} exhausted retries. Trying fallback model...`);
          break;
        }

        if (error instanceof AnalysisError && !isRateLimit && !isOverloaded) {
          throw error;
        }

        console.error("Gemini API Error:", error);
        
        if (error.message?.includes("SAFETY")) {
          throw new AnalysisError("The input contains content that triggered safety filters. Please check your compound names or structures.", "SAFETY_TRIGGERED");
        }

        if (isPermissionDenied) {
          throw new AnalysisError("Gemini API access denied (Permission Denied). Please check that your Gemini API key is valid and has the Generative Language API enabled with unrestricted domain/referrer access.", "PERMISSION_DENIED");
        }

        if (isRateLimit) {
          throw new AnalysisError("The daily request quota for the Gemini API has been reached or you are sending requests too quickly. Please wait 60 seconds and try again.", "QUOTA_EXCEEDED");
        }

        if (error.message?.includes("network") || error.message?.includes("fetch")) {
          throw new AnalysisError("A network connection issue was detected. Please check your internet connection and verify that the API is reachable.", "CONNECTION_ERROR");
        }
        
        if (isOverloaded) {
          throw new AnalysisError("The AI engine is currently experiencing high volume and is temporarily overloaded. Please try your request again in a few moments.", "MODEL_OVERLOADED");
        }

        if (error.message?.includes("API key not valid") || error.message?.includes("API_KEY_INVALID")) {
          throw new AnalysisError("The configured Gemini API key is invalid or has been revoked. Please verify your credentials in the Settings menu.", "CONFIG_ERROR");
        }

        // Default error with more context
        const errorMessage = error.message || "An analytical failure occurred while processing the molecular structures.";
        throw new AnalysisError(`${errorMessage}`, "UNKNOWN_ERROR");
      }
    }
  }

  throw new AnalysisError("Failed after trying available AI models. Please try again in a few moments.", "UNKNOWN_ERROR");
}
