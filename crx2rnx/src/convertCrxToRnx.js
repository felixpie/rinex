const fs = require('fs');
const path = require('path');
const wasm = require('../pkg');

async function convertCrxToRnx(inputFilePath, outputFilePath) {
  // Read the input compact RINEX file
  const inputFile = fs.readFileSync(inputFilePath);

  // Convert the compact RINEX to RINEX format using WebAssembly
  const rinexData = wasm.convert_crx_to_rnx(inputFile);

  // Write the output RINEX file
  fs.writeFileSync(outputFilePath, rinexData);
}

module.exports = convertCrxToRnx;
