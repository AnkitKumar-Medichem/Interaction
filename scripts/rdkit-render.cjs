const fs = require('fs');
const path = require('path');

let rdkitPromise = null;

function getRDKit() {
  if (!rdkitPromise) {
    const rdkitCode = fs.readFileSync(path.resolve(__dirname, '../public/RDKit_minimal.js'), 'utf8');
    const fn = new Function('require', '__dirname', '__filename', rdkitCode + '; return initRDKitModule;');
    const init = fn(require, __dirname, __filename);
    rdkitPromise = init({
      locateFile: () => path.resolve(__dirname, '../public/RDKit_minimal.wasm')
    });
  }
  return rdkitPromise;
}

async function main() {
  const smiles = process.argv[2];
  if (!smiles) {
    console.log(JSON.stringify({ error: 'No SMILES provided' }));
    return;
  }

  const width = parseInt(process.argv[3], 10) || 240;
  const height = parseInt(process.argv[4], 10) || 240;

  try {
    const rdkit = await getRDKit();
    const mol = rdkit.get_mol(smiles);
    if (!mol || !mol.is_valid()) {
      console.log(JSON.stringify({ error: 'Invalid SMILES', valid: false }));
      return;
    }

    const svg = mol.get_svg(width, height);
    let descriptors = {};
    try {
      const descStr = mol.get_descriptors();
      descriptors = JSON.parse(descStr);
    } catch (e) {
      descriptors = {};
    }

    mol.delete();

    console.log(JSON.stringify({
      valid: true,
      svg,
      descriptors: {
        MolWt: descriptors.exactmw || descriptors.amw || 0,
        MolLogP: descriptors.CrippenClogP || 0,
        TPSA: descriptors.tpsa || 0,
        NumRotatableBonds: descriptors.NumRotatableBonds || 0
      }
    }));
  } catch (err) {
    console.log(JSON.stringify({ error: err.message, valid: false }));
  }
}

main();
