import "dotenv/config";
import express, { Request, Response } from "express";
import path from "path";
import { GoogleGenAI, Type } from "@google/genai";
import { parse as parsePartial } from "partial-json";
import { createServer as createViteServer } from "vite";
import { generateComputationalPrediction, lookupCompoundSmiles } from "./src/lib/reaction-engine";
import { generateInteractionHeatmap } from "./src/lib/heatmap-engine";
import { renderSeabornHeatmap } from "./src/lib/python-seaborn";

const app = express();
const PORT = 3000;

app.use(express.json({ limit: "1mb" }));

// In-memory sliding-window IP rate limiter to prevent abuse and quota exhaustion
interface RateRecord {
  count: number;
  resetTime: number;
}
const rateLimits = new Map<string, RateRecord>();

function checkRateLimit(ip: string, limit: number, windowMs: number): boolean {
  const now = Date.now();
  const record = rateLimits.get(ip);
  if (!record || now > record.resetTime) {
    rateLimits.set(ip, { count: 1, resetTime: now + windowMs });
    return true;
  }
  if (record.count >= limit) {
    return false;
  }
  record.count++;
  return true;
}

// Clean up expired rate records every 5 minutes
setInterval(() => {
  const now = Date.now();
  for (const [key, record] of rateLimits.entries()) {
    if (now > record.resetTime) {
      rateLimits.delete(key);
    }
  }
}, 300000);

// Helper to classify Gemini API errors and determine if they are fatal or retryable
function classifyGeminiError(err: any): { type: string; message: string; isFatal: boolean } {
  const errMsg = err?.message || (typeof err === "string" ? err : JSON.stringify(err));
  const lower = errMsg.toLowerCase();

  if (
    errMsg.includes("API_KEY_INVALID") ||
    errMsg.includes("API key not valid") ||
    lower.includes("api key not valid") ||
    lower.includes("invalid api key")
  ) {
    return {
      type: "CONFIG_ERROR",
      message: "Your Gemini API key is invalid. The key currently saved in environment secrets appears to be an invalid key or a Firebase Web API key instead of a Gemini API key. Please generate a valid Gemini API key at https://aistudio.google.com/apikey and update the GEMINI_API_KEY secret in the AI Studio Settings panel.",
      isFatal: true
    };
  }

  if (errMsg.includes("PERMISSION_DENIED") || lower.includes("permission denied") || errMsg.includes("403")) {
    return {
      type: "PERMISSION_DENIED",
      message: "Gemini API access denied (Permission Denied). Please verify that your Gemini API key has the Generative Language API enabled with unrestricted domain access.",
      isFatal: true
    };
  }

  if (errMsg.includes("RESOURCE_EXHAUSTED") || errMsg.includes("429") || lower.includes("quota") || lower.includes("rate limit")) {
    return {
      type: "QUOTA_EXCEEDED",
      message: "The Gemini API request quota has been reached or requests are being sent too quickly. Please wait 60 seconds and try again.",
      isFatal: false
    };
  }

  if (errMsg.includes("SAFETY") || lower.includes("safety filter")) {
    return {
      type: "SAFETY_TRIGGERED",
      message: "The input chemical query triggered safety filters. Please verify that input compound names and SMILES do not include prohibited materials.",
      isFatal: true
    };
  }

  if (errMsg.includes("UNAVAILABLE") || errMsg.includes("503") || lower.includes("overloaded") || lower.includes("high demand")) {
    return {
      type: "MODEL_OVERLOADED",
      message: "The AI service is temporarily experiencing high volume and is overloaded. Please try again in a few moments.",
      isFatal: false
    };
  }

  return {
    type: "UNKNOWN_ERROR",
    message: errMsg || "An analytical failure occurred while generating chemical predictions.",
    isFatal: false
  };
}

// Track API key operational status in memory
let geminiKeyStatus: "untested" | "valid" | "invalid" = "untested";
let lastCheckedApiKey = "";
let probePromise: Promise<boolean> | null = null;

async function probeGeminiApiKeySilently(key: string): Promise<boolean> {
  if (!key || key === "MISSING_KEY" || key === "INVALID_OR_MISSING_KEY" || key.trim().length === 0) {
    return false;
  }
  try {
    const ai = new GoogleGenAI({ apiKey: key });
    await ai.models.generateContent({
      model: "gemini-2.5-flash",
      contents: "ping",
      config: { maxOutputTokens: 1 }
    });
    return true;
  } catch {
    return false;
  }
}

async function isGeminiKeyConfiguredAndValid(): Promise<boolean> {
  const apiKey = process.env.GEMINI_API_KEY;
  if (!apiKey || apiKey === "MISSING_KEY" || apiKey === "INVALID_OR_MISSING_KEY" || apiKey.trim().length === 0) {
    geminiKeyStatus = "invalid";
    return false;
  }
  if (apiKey !== lastCheckedApiKey || geminiKeyStatus === "untested") {
    lastCheckedApiKey = apiKey;
    if (!probePromise) {
      probePromise = probeGeminiApiKeySilently(apiKey).then((ok) => {
        geminiKeyStatus = ok ? "valid" : "invalid";
        probePromise = null;
        return ok;
      });
    }
    return await probePromise;
  }
  return geminiKeyStatus === "valid";
}

// Background validation probe on boot
if (process.env.GEMINI_API_KEY) {
  isGeminiKeyConfiguredAndValid().catch(() => {});
}

const getAiClient = () => {
  const apiKey = process.env.GEMINI_API_KEY;
  if (!apiKey || apiKey === "MISSING_KEY" || apiKey === "INVALID_OR_MISSING_KEY" || apiKey.trim().length === 0) {
    throw new Error("GEMINI_API_KEY is not configured on the server.");
  }
  return new GoogleGenAI({ apiKey });
};

// ==========================================
// 1. Prediction Endpoint (Streaming SSE)
// ==========================================
app.post("/api/predict", async (req: Request, res: Response) => {
  const clientIp = req.ip || req.socket.remoteAddress || "unknown";
  
  // Rate limit: Max 20 predictions per minute per IP
  if (!checkRateLimit(`predict:${clientIp}`, 20, 60000)) {
    return res.status(429).json({ error: "Rate limit exceeded. Please wait 60 seconds before making more predictions." });
  }

  const { inputs, method = "Heuristic" } = req.body;

  // Strict input validation
  if (!Array.isArray(inputs) || inputs.length === 0 || inputs.length > 5) {
    return res.status(400).json({ error: "Invalid inputs: Provide between 1 and 5 compounds." });
  }

  for (const input of inputs) {
    if (!input || typeof input.value !== "string" || input.value.trim().length === 0 || input.value.length > 1000) {
      return res.status(400).json({ error: "Invalid compound value. Must be a non-empty string under 1000 characters." });
    }
    if (input.type !== "Name" && input.type !== "SMILES") {
      return res.status(400).json({ error: "Invalid input type. Must be 'Name' or 'SMILES'." });
    }
  }

  const validMethods = ["Boltzmann", "Heuristic", "Both"];
  if (!validMethods.includes(method)) {
    return res.status(400).json({ error: "Invalid prediction method." });
  }

  // Setup Server-Sent Events (SSE) for secure real-time streaming to the browser
  res.setHeader("Content-Type", "text/event-stream");
  res.setHeader("Cache-Control", "no-cache, no-transform");
  res.setHeader("Connection", "keep-alive");
  res.flushHeaders();

  const sendSse = (event: string, data: any) => {
    res.write(`event: ${event}\ndata: ${JSON.stringify(data)}\n\n`);
  };

  try {
    const compoundsInfo = inputs
      .map((input: any, i: number) => {
        let desc = "";
        if (input.descriptors) {
          desc = ` [Calculated Specs: MolWt: ${input.descriptors.MolWt?.toFixed(2) || "N/A"}, LogP: ${input.descriptors.MolLogP?.toFixed(2) || "N/A"}, TPSA: ${input.descriptors.TPSA?.toFixed(2) || "N/A"}, Rotatable Bonds: ${input.descriptors.NumRotatableBonds ?? "N/A"}]`;
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
      origin: { type: Type.STRING, description: "Which specific compound(s) this product originated from, such as 'Compound 1' or 'Compound 1 and Compound 2'." },
      probability: { 
        type: Type.NUMBER, 
        description: "Primary probability of formation as a decimal between 0.0 and 1.0. Used for ranking." 
      },
      condition: { 
        type: Type.STRING, 
        enum: ["Oxidation", "Acidic Hydrolysis", "Basic Hydrolysis", "Hydrolysis", "Photodegradation", "Thermal Degradation"] 
      },
      source: { 
        type: Type.STRING, 
        enum: ["Stress degradation", "Interaction with other compound"] 
      },
      mechanismExplanation: {
        type: Type.STRING,
        description: "Brief explanation of the interaction mechanism including pH change, oxidation, complexation, adsorption, or precipitation. Mention effects on stability or release kinetics."
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

    let fullText = "";
    let lastClassifiedError: { type: string; message: string; isFatal: boolean } | null = null;

    const useAi = await isGeminiKeyConfiguredAndValid();

    if (useAi) {
      const modelsToTry = ["gemini-3.8-flash", "gemini-flash-latest", "gemini-3.1-flash-lite"];
      const maxRetries = 1;

      modelLoop: for (let mIdx = 0; mIdx < modelsToTry.length; mIdx++) {
        const activeModel = modelsToTry[mIdx];
        let attempt = 0;
        let success = false;

        while (attempt <= maxRetries && !success) {
          try {
            const ai = getAiClient();
            const responseStream = await ai.models.generateContentStream({
              model: activeModel,
              contents: `Predict and evaluate the chemical interaction and reaction products of Compound 1 in the following mixture using the ${method === "Both" ? "Heuristic AND Boltzmann" : method}-based approach:\n${compoundsInfo}`,
              config: {
                temperature: 0.1,
                systemInstruction: `You are an expert computational chemist, cheminformatician, and reaction mechanism evaluator.
CRITICAL CALCULATION RULE: All calculations and predictions MUST be carried out based on:
1. Systematic identification of every functional group in the input molecule (Compound 1) and any secondary compounds (Compounds 2-5).
2. The specific chemical reactivity of those identified functional groups against:
   - Acidic stress (including A_Ac2 solvolysis, hydronium protonation, acid-catalyzed dehydration)
   - Basic stress (including B_Ac2 saponification, hydroxide nucleophilic attack, deprotonation)
   - Hydrolysis (neutral water solvolysis across labile linkages under humidity)
   - Photolytic stress (UV chromophore excitation, photo-Fries rearrangement, Norrish cleavage, photo-oxidation)
   - Thermal stress (pyrolysis, decarboxylation, syn-elimination, thermal condensation)
   - Oxidative stress (single-electron transfer, phenoxy/anilinyl radicals, S-oxidation to sulfoxide/sulfone, N-oxidation)
   - AND/OR cross-reactions with functional groups of secondary compound(s) (including transamidation, Maillard Schiff base with reducing sugars, transesterification, chelation/salt formation).

Analytical Framework${method === "Both" ? "s" : ""}:
${method === "Heuristic" || method === "Both" ? "1. HEURISTIC ANALYSIS: Based on expert chemical reasoning, functional group reactive sites, and known reaction kinetics." : ""}
${method === "Boltzmann" || method === "Both" ? `${method === "Both" ? "2." : "1."} BOLTZMANN ANALYSIS: Based on thermodynamic stability and calculated relative formation free energy (ΔG in kcal/mol at 298.15K) via Boltzmann distribution.` : ""}

${method === "Both" ? "When 'Both' is selected, perform both analyses independently for each predicted product to provide a comparative perspective." : ""}

First, identify the chemical structures and functional groups correctly for ALL provided compounds.
For each compound, provide:
- Identified name (Echo the exact user-specified name or standard chemical name).
- SMILES string (MUST be valid, canonical SMILES obeying valency rules).
- List of key functional group features.
- List of specific reactive interaction sites.

Predict ONLY the TOP 5 most significant reaction byproducts, degradation products, or interaction adducts derived from Compound 1.
For each product, you MUST specify:
- Whether it forms from 'Direct degradation' or 'Interaction with other compound'.
- Which specific condition it forms under ('Acidic Hydrolysis', 'Basic Hydrolysis', 'Hydrolysis', 'Photodegradation', 'Thermal Degradation', 'Oxidation').
- IUPAC name of the NEW PRODUCT.
- SMILES string of the new product (valid canonical SMILES obeying valences).
- A detailed explanation of the underlying mechanism explicitly citing the reacting functional group(s) (mechanismExplanation).
- ${probabilityInstruction}

IMPORTANT: Probabilities MUST be realistic estimates between 0.01 and 0.99.
Rank the products by their calculated probability descending. Do not return more than 5 products.`,
                responseMimeType: "application/json",
                responseSchema: {
                  type: Type.OBJECT,
                  properties: {
                    chainOfThought: {
                      type: Type.STRING,
                      description: `Step-by-step chemical reasoning: detail the functional group identification of Compound 1, evaluate condition-by-condition reactivity (acidic, basic, hydrolysis, photolytic, thermal, oxidative), calculate cross-reactivity with secondary compounds, and estimate formation energies/probabilities.`
                    },
                    compounds: {
                      type: Type.ARRAY,
                      items: {
                        type: Type.OBJECT,
                        properties: {
                          name: { type: Type.STRING },
                          smiles: { type: Type.STRING },
                          features: { type: Type.ARRAY, items: { type: Type.STRING } },
                          interactionSites: { type: Type.ARRAY, items: { type: Type.STRING } }
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

            fullText = "";
            for await (const chunk of responseStream) {
              if (chunk.text) {
                fullText += chunk.text;
                try {
                  const partial = parsePartial(fullText);
                  if (partial) {
                    sendSse("chunk", partial);
                  }
                } catch (_) {}
              }
            }

            if (fullText) {
              const parsed = JSON.parse(fullText);
              if (parsed && Array.isArray(parsed.degradationImpurities)) {
                parsed.degradationImpurities = parsed.degradationImpurities.slice(0, 5);
              }
              sendSse("complete", parsed);
              success = true;
              geminiKeyStatus = "valid";
              break modelLoop;
            }
          } catch (err: any) {
            const classified = classifyGeminiError(err);
            lastClassifiedError = classified;

            if (classified.isFatal || classified.type === "CONFIG_ERROR" || classified.type === "PERMISSION_DENIED") {
              geminiKeyStatus = "invalid";
              break modelLoop;
            }

            attempt++;
            if (attempt <= maxRetries) {
              await new Promise(r => setTimeout(r, 1000 * attempt));
            }
          }
        }
      }
    }

    if (!fullText) {
      const engineResult = generateComputationalPrediction(inputs, method);
      
      // Send progressive thinking chunk
      sendSse("chunk", {
        chainOfThought: engineResult.chainOfThought,
        compounds: engineResult.compounds
      });

      // Complete report
      sendSse("complete", engineResult);
      return;
    }
  } catch (error: any) {
    try {
      const fallbackResult = generateComputationalPrediction(req.body.inputs || [], req.body.method || "Both");
      sendSse("complete", fallbackResult);
    } catch (engineErr) {
      const classified = classifyGeminiError(error);
      sendSse("error", { type: classified.type, message: classified.message });
    }
  } finally {
    res.end();
  }
});

// ==========================================
// 2. SMILES Remediation Endpoint
// ==========================================
app.post("/api/remediate-smiles", async (req: Request, res: Response) => {
  const clientIp = req.ip || req.socket.remoteAddress || "unknown";
  
  // Rate limit: Max 60 remediations per minute per IP
  if (!checkRateLimit(`remediate:${clientIp}`, 60, 60000)) {
    return res.status(429).json({ error: "Rate limit exceeded for SMILES lookup." });
  }

  const { name } = req.body;
  if (!name || typeof name !== "string" || name.trim().length === 0 || name.length > 200) {
    return res.status(400).json({ error: "Invalid compound name." });
  }

  // 1. High-speed lookup in local pharmaceutical repository
  const localSmiles = lookupCompoundSmiles(name);
  if (localSmiles) {
    return res.json({ smiles: localSmiles });
  }

  // 2. Check if external AI key is configured and valid
  const hasValidKey = await isGeminiKeyConfiguredAndValid();
  if (!hasValidKey) {
    return res.json({ smiles: null });
  }

  // 3. AI model lookup fallback
  try {
    const ai = getAiClient();
    const response = await ai.models.generateContent({
      model: "gemini-3.8-flash",
      contents: `Provide the valid, canonical SMILES string for the compound named "${name.trim()}". Return ONLY the SMILES string.`,
      config: {
        responseMimeType: "application/json",
        responseSchema: {
          type: Type.OBJECT,
          properties: {
            smiles: { type: Type.STRING }
          },
          required: ["smiles"]
        }
      }
    });

    const result = JSON.parse(response.text || "{}");
    return res.json({ smiles: result.smiles || null });
  } catch (error: any) {
    const classified = classifyGeminiError(error);
    if (classified.isFatal || classified.type === "CONFIG_ERROR" || classified.type === "PERMISSION_DENIED") {
      geminiKeyStatus = "invalid";
    }
    return res.json({ smiles: null });
  }
});

// ==========================================
// 3. Reactive Centers Interaction Heatmap
// ==========================================
app.post("/api/interaction-heatmap", async (req: Request, res: Response) => {
  try {
    const { compounds, cmap = "warmcool", title, conditions } = req.body;
    if (!Array.isArray(compounds) || compounds.length === 0) {
      return res.status(400).json({ error: "At least one compound is required for interaction heatmap." });
    }

    const result = generateInteractionHeatmap({ compounds, cmap, title, conditions });

    // Render publication-grade Seaborn heatmap using Python backend
    try {
      const chartTitle =
        title ||
        (result.compound2Name !== "Intramolecular"
          ? `Stress Degradation & Incompatibility Heatmap: ${result.compound1Name} & ${result.compound2Name}`
          : `Stress Degradation Heatmap: ${result.compound1Name}`);

      const seabornPng = await renderSeabornHeatmap({
        matrix: result.matrix,
        rowLabels: result.rowLabels,
        colLabels: result.colLabels,
        title: chartTitle,
        cmap,
      });

      if (seabornPng) {
        result.heatmapBase64 = seabornPng;
      }
    } catch (seabornErr) {
      console.warn("Seaborn heatmap generation fallback to SVG:", seabornErr);
    }

    return res.json(result);
  } catch (err: any) {
    return res.status(500).json({ error: err?.message || "Internal server error" });
  }
});

// ==========================================
// 4. Vite Middleware & Asset Serving
// ==========================================
async function startServer() {
  if (process.env.NODE_ENV !== "production") {
    const vite = await createViteServer({
      server: { middlewareMode: true },
      appType: "spa",
    });
    app.use(vite.middlewares);
  } else {
    const distPath = path.join(process.cwd(), "dist");
    app.use(express.static(distPath));
    app.get("*", (_req, res) => {
      res.sendFile(path.join(distPath, "index.html"));
    });
  }

  app.listen(PORT, "0.0.0.0", () => {
    console.log(`Server running securely on port ${PORT}`);
  });
}

startServer();
