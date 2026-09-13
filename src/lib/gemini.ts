import { validateSmiles, MolecularDescriptors } from "./rdkit";

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

/**
 * Securely proxies prediction requests to the server backend.
 * Protects the private GEMINI_API_KEY from exposure to browser inspection or network traces.
 */
export async function predictInteraction(
  inputs: CompoundInput[],
  method: PredictionMethod = "Heuristic",
  onChunk?: (partialResult: Partial<PredictionResult>) => void
): Promise<PredictionResult> {
  // Pre-validation with client RDKit
  for (const input of inputs) {
    if (input.type === "SMILES") {
      const smiles = input.value.trim();
      const validation = await validateSmiles(smiles);
      if (!validation.isValid) {
        throw new AnalysisError(validation.error || "Invalid chemical structure", "INVALID_SMILES");
      }
      if (validation.canonicalSmiles) {
        input.value = validation.canonicalSmiles;
      }
    }
  }

  const response = await fetch("/api/predict", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "Accept": "text/event-stream"
    },
    body: JSON.stringify({ inputs, method })
  });

  if (!response.ok) {
    let errorMsg = `Server error (${response.status})`;
    let errorType = "UNKNOWN_ERROR";
    try {
      const errJson = await response.json();
      if (errJson.error) {
        errorMsg = errJson.error;
        if (response.status === 429) errorType = "QUOTA_EXCEEDED";
        else if (response.status === 400) errorType = "INVALID_INPUT";
      }
    } catch (_) {}
    throw new AnalysisError(errorMsg, errorType);
  }

  const reader = response.body?.getReader();
  if (!reader) {
    throw new AnalysisError("Unable to establish streaming response connection.", "CONNECTION_ERROR");
  }

  const decoder = new TextDecoder();
  let buffer = "";
  let finalResult: PredictionResult | null = null;

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    buffer += decoder.decode(value, { stream: true });
    const lines = buffer.split("\n");
    buffer = lines.pop() || "";

    let currentEvent = "message";
    for (const line of lines) {
      if (line.startsWith("event: ")) {
        currentEvent = line.slice(7).trim();
      } else if (line.startsWith("data: ")) {
        const dataStr = line.slice(6).trim();
        if (!dataStr) continue;
        try {
          const payload = JSON.parse(dataStr);
          if (currentEvent === "chunk" && onChunk) {
            onChunk(payload as Partial<PredictionResult>);
          } else if (currentEvent === "complete") {
            finalResult = payload as PredictionResult;
          } else if (currentEvent === "error") {
            throw new AnalysisError(payload.message || "An analytical failure occurred.", payload.type || "SERVER_ERROR");
          }
        } catch (e) {
          if (e instanceof AnalysisError) throw e;
        }
      }
    }
  }

  if (!finalResult) {
    throw new AnalysisError("Analysis completed without a complete response report. Please try again.", "EMPTY_RESPONSE");
  }

  return finalResult;
}

