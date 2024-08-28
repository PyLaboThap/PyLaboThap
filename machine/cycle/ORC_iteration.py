# class RC(Cycle):
#     def __init__(self, fluid=None):
#         super().__init__(fluid)  # Ensure Cycle's initialization is done

#     def solve(self, pressure_guesses, residuals, tolerance=1e-3, max_iterations=1000):
#         def objective(x):
#             P_ev, rp = x
#             P_cd = rp * P_ev
#             print(x, 'x')
#             # Step 1: Set the initial guesses and calculate related properties
#             self.set_cycle_guesses(P_ev, P_cd)

#             # Access the evaporator and condenser components directly
#             evaporator = self.get_component("Evaporator")
#             condenser = self.get_component("Condenser")
#             pump = self.get_component("Pump")
#             expander = self.get_component("Expander")

#             # Step 2: Evaluate initial enthalpies from the component ports
#             print(evaporator)
#             # print(pump.model.su.p)
#             # print(pump.model.su.h)
#             print(condenser.model.ex_wf.p)
#             print(evaporator.model.ex_wf.h)
#             h_ex_ev_init = evaporator.model.ex_wf.h
#             h_ex_cd_init = condenser.model.ex_wf.h

#             print(self.components.values())
#             # Step 3: Solve the cycle with the guesses
#             for component in self.components.values():
#                 component.solve()
#                 print(component)

#             print(condenser.model.ex_wf.p)
#             print(pump.model.su.p)
#             print(pump.model.ex.p)
#             # Step 4: Recalculate the enthalpies after solving
#             h_ex_ev_final = evaporator.model.ex_wf.h
#             h_ex_cd_final = condenser.model.ex_wf.h

#             # Step 5: Compute the residuals as differences in enthalpies
#             residual_h_ex_ev = h_ex_ev_init - h_ex_ev_final
#             residual_h_ex_cd = h_ex_cd_init - h_ex_cd_final
#             print(residual_h_ex_ev, 'residual_h_ex_ev')
#             print(residual_h_ex_cd, 'residual_h_ex_cd')
#             res = np.sqrt(residual_h_ex_ev ** 2 + residual_h_ex_cd ** 2)
#             print(res, 'res')
#             # res = [residual_h_ex_ev, residual_h_ex_cd]
#             return res

#         # Set up initial guesses for P_ev and P_cd
#         x0 = [pressure_guesses["P_ev"], 3]
#         bounds = [(1e4, 5e6), (1, 10)]  # Example: (10 kPa to 5 MPa)

#         # Optimize using scipy's minimize function
#         result = minimize(objective, x0, bounds=bounds, tol=tolerance, options={'maxiter': max_iterations})
#         # result = fsolve(objective, x0)
#         print(result, 'result')

#         if result.success:
#             optimized_pressures = result.x
#             print(f"Optimization succeeded: P_ev = {optimized_pressures[0]}, P_cd = {optimized_pressures[1]}")
#         else:
#             print("Optimization failed.")

#     def set_cycle_guesses(self, P_ev_guess, P_cd_guess):
#         """
#         Set the initial guesses for the evaporator and condenser pressures and calculate related properties.
#         """
#         subcooling = 5  # Adjust as needed
#         superheating = 5  # Adjust as needed

#         # Step 1: Calculate temperature at the pump inlet based on P_cd_guess and subcooling
#         T_su_pp = PropsSI('T', 'P', P_cd_guess, 'Q', 0, self.fluid) - subcooling
#         self.set_guess(target="Pump:su", T=T_su_pp, P=P_cd_guess, fluid=self.fluid)
#         # self.set_guess(target="Pump:su", P=P_cd_guess)

#         # Step 2: Set guessed pressure at the expander inlet (evaporation pressure)
#         self.set_guess(target="Pump:ex", P=P_ev_guess)
#         T_ex_ev = PropsSI('T', 'P', P_ev_guess, 'Q', 1, self.fluid) + superheating
#         self.set_guess(target="Evaporator:ex_wf", P=P_ev_guess, T=T_ex_ev)
#         self.set_guess(target="Expander:ex", P=P_cd_guess)