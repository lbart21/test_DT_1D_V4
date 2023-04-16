model = "CEAGas"

CEAGas = {
  mixtureName = 'h2-o2-n2-9sp',
  speciesList = {"H2", "O2", "N2", "O", "H", "OH", "H2O2", "H2O", "HO2"},
  reactants = {H2=0.5, O2=0.5, N2=0.0, O=0.0, H=0.0, OH=0.0, H2O=0.0, HO2=0.0, H2O2=0.0},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}