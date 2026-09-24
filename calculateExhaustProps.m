function [MW_exhaust, R_spec] = calculateExhaustProps(OF_ratio, w_ethanol, w_water)
% CALCULATEEXHAUSTPROPS Computes mean molecular weight and gas constant.
% Inputs:
%   OF_ratio  - Mass ratio of N2O to fuel blend (e.g., 2.1:1)
%   w_ethanol - Weight fraction of ethanol in fuel (e.g., 0.70)
%   w_water   - Weight fraction of water in fuel (e.g., 0.30)
% Outputs:
%   MW_exhaust - Mean molecular weight of exhaust [kg/mol]
%   R_spec     - Specific gas constant [J/(kg*K)]
   % Molar masses of reactants [kg/mol]
   MW_N2O     = 0.044013;
   MW_Ethanol = 0.046069;
   MW_Water   = 0.018015;
   % Mass balance based on 1.0 kg of fuel
   m_fuel  = 1.0;
   m_N2O   = OF_ratio * m_fuel;
   m_total = m_fuel + m_N2O;
   % Moles of reactants
   n_ethanol = (m_fuel * w_ethanol) / MW_Ethanol;
   n_water   = (m_fuel * w_water)   / MW_Water;
   n_N2O     = m_N2O / MW_N2O;
   % Elemental atom balance
   n_C = 2 * n_ethanol;
   n_H = 6 * n_ethanol + 2 * n_water;
   n_O = 1 * n_ethanol + 1 * n_water + 1 * n_N2O;
   n_N = 2 * n_N2O;
   % Equilibrium product estimation
   n_N2  = n_N / 2;
   n_CO2 = 0.12 * n_C;
   n_CO  = n_C - n_CO2;
   n_H2O = n_O - (2 * n_CO2 + n_CO);
   n_H2  = (n_H - 2 * n_H2O) / 2;
   n_total_gas = n_N2 + n_CO + n_CO2 + n_H2O + n_H2;
   % Output calculations
   MW_exhaust = m_total / n_total_gas;         % [kg/mol]
   Ru         = 8.3144626185324;              % Universal gas constant [J/(mol*K)]
   R_spec     = Ru / MW_exhaust;             % Specific gas constant [J/(kg*K)]
end

