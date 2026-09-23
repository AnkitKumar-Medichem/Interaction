import { MolecularDescriptors } from "./lib/rdkit";

export interface CompoundInfo {
  name: string;
  smiles: string;
  features: string[];
  interactionSites: string[];
  molecularDescriptors?: MolecularDescriptors;
}

export interface FunctionalGroupReactivity {
  groupName: string;
  category: string;
  smilesFragment: string;
  reactiveSite: string;
  acidic: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "Resistant"; mechanism: string };
  basic: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "Resistant"; mechanism: string };
  hydrolysis: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "Resistant"; mechanism: string };
  photolytic: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "Resistant"; mechanism: string };
  thermal: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "Resistant"; mechanism: string };
  oxidative: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "Resistant"; mechanism: string };
  secondaryInteraction?: { vulnerability: "Critical" | "High" | "Moderate" | "Low" | "None" | "Resistant"; partnerGroup?: string; mechanism: string };
}

export interface PredictionResult {
  chainOfThought: string;
  compounds: CompoundInfo[];
  interactionType: "Physical" | "Chemical" | "None";
  mechanism: string;
  functionalGroupAnalysis?: FunctionalGroupReactivity[];
  degradationImpurities: {
    iupacName: string;
    smiles: string;
    structureDescription: string;
    origin: string;
    probability: number; // Primary probability for ranking
    probabilityHeuristic?: number;
    probabilityBoltzmann?: number;
    relativeEnergy?: number;
    condition: "Oxidation" | "Acidic Hydrolysis" | "Basic Hydrolysis" | "Hydrolysis" | "Photodegradation" | "Thermal Degradation" | string;
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
