import { GoogleGenAI, Type } from "@google/genai";
import { updateCompoundSmiles } from './firestore-utils';

const getApiKey = () => {
  if (typeof window !== "undefined") {
    const win = window as any;
    if (win.__GEMINI_API_KEY__ && win.__GEMINI_API_KEY__ !== "MISSING_KEY") {
      return win.__GEMINI_API_KEY__;
    }
    const params = new URLSearchParams(window.location.search);
    const urlKey = params.get("gemini_key") || params.get("api_key");
    if (urlKey) return urlKey;
    try {
      const storedKey = localStorage.getItem("GEMINI_API_KEY");
      if (storedKey && storedKey !== "MISSING_KEY") return storedKey;
    } catch (_) {}
  }
  const key = process.env.GEMINI_API_KEY || process.env.VITE_GEMINI_API_KEY;
  return key && key !== "MISSING_KEY" ? key : null;
};

const getAiClient = () => {
  return new GoogleGenAI({ apiKey: getApiKey() || "" });
};

/**
 * Dynamically resolves a missing SMILES string using Gemini and updates the Firestore database.
 */
export async function remediateCompoundSmiles(docId: string, name: string): Promise<string | null> {
  const apiKey = getApiKey();
  if (!apiKey) {
    console.error("Gemini API key is required to resolve missing SMILES.");
    return null;
  }

  const maxRetries = 2;
  let attempt = 0;

  while (attempt <= maxRetries) {
    try {
      const ai = getAiClient();
      const response = await ai.models.generateContent({
        model: "gemini-3.1-flash-lite",
        contents: `Provide the valid, canonical SMILES string for the compound named "${name}". Return ONLY the SMILES string.`,
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
      const smiles = result.smiles;
      
      if (smiles) {
        await updateCompoundSmiles(docId, smiles);
        return smiles;
      }
      return null;
    } catch (error: any) {
      const isRateLimit = error.message?.includes("quota") || 
                          error.message?.includes("429") || 
                          error.status === 429 ||
                          error.message?.toLowerCase().includes("rate limit");
                          
      if (attempt < maxRetries && isRateLimit) {
        attempt++;
        const delayMs = attempt * 2000 + Math.random() * 1000;
        console.warn(`Rate limit fetching SMILES for ${name}. Retrying in ${Math.round(delayMs/1000)}s...`);
        await new Promise(resolve => setTimeout(resolve, delayMs));
        continue;
      }
      
      console.error(`Failed to remediate SMILES for ${name}:`, error);
      return null;
    }
  }
  return null;
}
