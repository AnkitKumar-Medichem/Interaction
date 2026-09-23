import { validateSmiles } from "./rdkit";
import { PredictionResult, PredictionMethod, CompoundInput, AnalysisError } from "../types";

export * from "../types";

/**
 * Predict chemical interaction products, stress degradants, and reaction mechanisms
 * via the local deterministic computational reaction engine.
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
  let currentEvent = "message";
  let lastPartial: Partial<PredictionResult> | null = null;

  const processLine = (line: string) => {
    const trimmed = line.trim();
    if (!trimmed) return;

    if (trimmed.startsWith("event: ")) {
      currentEvent = trimmed.slice(7).trim();
    } else if (trimmed.startsWith("data: ")) {
      const dataStr = trimmed.slice(6).trim();
      if (!dataStr) return;
      try {
        const payload = JSON.parse(dataStr);
        if (currentEvent === "chunk") {
          lastPartial = payload as Partial<PredictionResult>;
          if (onChunk) {
            onChunk(payload as Partial<PredictionResult>);
          }
        } else if (currentEvent === "complete") {
          finalResult = payload as PredictionResult;
        } else if (currentEvent === "error") {
          throw new AnalysisError(payload.message || "An analytical failure occurred.", payload.type || "SERVER_ERROR");
        }
      } catch (e) {
        if (e instanceof AnalysisError) throw e;
      }
    }
  };

  while (true) {
    const { done, value } = await reader.read();
    if (done) break;

    buffer += decoder.decode(value, { stream: true });
    const lines = buffer.split("\n");
    buffer = lines.pop() ?? "";

    for (const line of lines) {
      processLine(line);
    }
  }

  buffer += decoder.decode();
  if (buffer.trim()) {
    const remainingLines = buffer.split("\n");
    for (const line of remainingLines) {
      processLine(line);
    }
  }

  if (!finalResult && lastPartial && lastPartial.compounds && lastPartial.mechanism) {
    finalResult = {
      chainOfThought: lastPartial.chainOfThought || "",
      compounds: lastPartial.compounds || [],
      interactionType: lastPartial.interactionType || "None",
      mechanism: lastPartial.mechanism || "",
      degradationImpurities: lastPartial.degradationImpurities || []
    };
  }

  if (!finalResult) {
    throw new AnalysisError("Analysis completed without a complete response report. Please try again.", "EMPTY_RESPONSE");
  }

  return finalResult;
}
