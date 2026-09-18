import React, { useEffect, useState } from 'react';
import { getMoleculeSvg, sanitizeSmiles } from '../lib/rdkit';

interface ChemicalStructureProps {
  smiles: string;
  width?: number;
  height?: number;
  className?: string;
}

export const ChemicalStructure: React.FC<ChemicalStructureProps> = ({ 
  smiles, 
  width = 200, 
  height = 200,
  className = ""
}) => {
  const [svgContent, setSvgContent] = useState<string | null>(null);
  const [fallbackImgUrl, setFallbackImgUrl] = useState<string | null>(null);
  const [imgError, setImgError] = useState<boolean>(false);
  const [error, setError] = useState<boolean>(false);

  useEffect(() => {
    let mounted = true;
    const clean = sanitizeSmiles(smiles);
    
    if (!clean) {
      setError(true);
      setSvgContent(null);
      setFallbackImgUrl(null);
      return;
    }

    setError(false);
    setImgError(false);
    
    getMoleculeSvg(clean, width, height).then((svg) => {
      if (!mounted) return;
      if (svg) {
        setSvgContent(svg);
        setFallbackImgUrl(null);
        setError(false);
      } else {
        // Fall back to chemical structure depiction service (Cactus / PubChem)
        setSvgContent(null);
        setFallbackImgUrl(`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(clean)}/image`);
      }
    }).catch(() => {
      if (!mounted) return;
      setSvgContent(null);
      setFallbackImgUrl(`https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(clean)}/image`);
    });

    return () => { mounted = false; };
  }, [smiles, width, height]);

  const handleImgError = () => {
    const clean = sanitizeSmiles(smiles);
    if (!imgError && fallbackImgUrl && fallbackImgUrl.includes("cactus.nci.nih.gov")) {
      // Try PubChem fallback
      setImgError(true);
      setFallbackImgUrl(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(clean)}/PNG?record_type=2d&image_size=300x300`);
    } else {
      setError(true);
      setFallbackImgUrl(null);
    }
  };

  return (
    <div 
      className={`relative flex items-center justify-center bg-white rounded-lg overflow-hidden ${className}`}
      style={{ width, height }}
    >
      {error && (
        <div className="absolute inset-0 flex items-center justify-center bg-slate-50 text-[10px] text-slate-400 p-2 text-center">
          Structure unavailable
        </div>
      )}
      {!error && svgContent && (
        <div 
          className="w-full h-full flex items-center justify-center [&>svg]:max-w-full [&>svg]:max-h-full [&>svg]:w-auto [&>svg]:h-auto"
          dangerouslySetInnerHTML={{ __html: svgContent }} 
        />
      )}
      {!error && !svgContent && fallbackImgUrl && (
        <img
          src={fallbackImgUrl}
          alt="Chemical Structure"
          className="max-w-full max-h-full object-contain p-2"
          onError={handleImgError}
          referrerPolicy="no-referrer"
        />
      )}
    </div>
  );
};
