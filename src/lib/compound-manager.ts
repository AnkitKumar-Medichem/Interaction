import { GoogleGenAI, Type } from "@google/genai";
import { updateCompoundSmiles } from './firestore-utils';

const getApiKey = () => {
  const key = process.env.GEMINI_API_KEY || process.env.VITE_GEMINI_API_KEY;
  return key && key !== "MISSING_KEY" ? key : null;
};

const ai = new GoogleGenAI({ apiKey: getApiKey() || "" });

/**
 * Dynamically resolves a missing SMILES string using Gemini and updates the Firestore database.
 */
export async function remediateCompoundSmiles(docId: string, name: string): Promise<string | null> {
  const apiKey = getApiKey();
  if (!apiKey) {
    console.error("Gemini API key is required to resolve missing SMILES.");
    return null;
  }

  try {
    const response = await ai.models.generateContent({
      model: "gemini-2.5-pro",
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
  } catch (error) {
    console.error(`Failed to remediate SMILES for ${name}:`, error);
  }
  return null;
}
