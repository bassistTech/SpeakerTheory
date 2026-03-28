function generated_code_1() {
C_ms = dprint("C_ms = ",V_as/(Math.pow(S_d, 2)*Math.pow(c, 2)*rho))
M_ms = dprint("M_ms = ",Math.pow(S_d, 2)*Math.pow(c, 2)*rho/(V_as*Math.pow(w_s, 2)))
R_ms = dprint("R_ms = ",Math.pow(S_d, 2)*Math.pow(c, 2)*rho/(Q_ms*V_as*w_s))
BL = dprint("BL = ",Math.sqrt(R_e)*S_d*c*Math.sqrt(rho)/(Math.sqrt(Q_es)*Math.sqrt(V_as)*Math.sqrt(w_s)))
P_ref = dprint("P_ref = ",2.00000000000000e-5)
R_ref = dprint("R_ref = ",1.00000000000000)
}
function generated_code_2(omega) {
C_eff_sealed = 1/(N_d*Math.pow(S_d, 2)*Math.pow(c, 2)*rho/V_box + 1/C_ms)
C_eff_ported = 1/(N_d*Math.pow(S_d, 2)*Math.pow(c, 2)*Math.pow(omega, 2)*rho/(V_box*(Math.pow(omega, 2) - Math.pow(w_port, 2))) + 1/C_ms)
if (Fport == 0) {
      C_eff = C_eff_sealed
      }
      else {
      C_eff = C_eff_ported
      }
      V = VinExc
x_squared_real = Math.pow(BL, 2)*Math.pow(C_eff, 2)/(Math.pow(BL, 4)*Math.pow(C_eff, 2)*Math.pow(omega, 2) - 2.0*Math.pow(BL, 2)*Math.pow(C_eff, 2)*L_e*M_ms*Math.pow(omega, 4) + 2*Math.pow(BL, 2)*Math.pow(C_eff, 2)*R_e*R_ms*Math.pow(omega, 2) + 2.0*Math.pow(BL, 2)*C_eff*L_e*Math.pow(omega, 2) + 1.0*Math.pow(C_eff, 2)*Math.pow(L_e, 2)*Math.pow(M_ms, 2)*Math.pow(omega, 6) + 1.0*Math.pow(C_eff, 2)*Math.pow(L_e, 2)*Math.pow(R_ms, 2)*Math.pow(omega, 4) + Math.pow(C_eff, 2)*Math.pow(M_ms, 2)*Math.pow(R_e, 2)*Math.pow(omega, 4) + Math.pow(C_eff, 2)*Math.pow(R_e, 2)*Math.pow(R_ms, 2)*Math.pow(omega, 2) - 2.0*C_eff*Math.pow(L_e, 2)*M_ms*Math.pow(omega, 4) - 2*C_eff*M_ms*Math.pow(R_e, 2)*Math.pow(omega, 2) + 1.0*Math.pow(L_e, 2)*Math.pow(omega, 2) + Math.pow(R_e, 2))
xabs = Math.sqrt(x_squared_real)
port_cone_frac = N_d*S_d*Math.pow(w_port, 2)/(S_port*(Math.pow(omega, 2) - Math.pow(w_port, 2)))
sound_pressure_ratio = (1/2)*S_d*Math.pow(omega, 2)*rho/(Math.PI*P_ref*R_ref)
}
