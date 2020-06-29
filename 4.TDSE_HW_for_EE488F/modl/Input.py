# Input parameter
Lx = 50            # Box Length (Angstrom = 10^-10 m)
ngx = 1001             # Number of grid points (spacing = L/(n-1))
nstep = 200            # Number of time steps * 0.1fs
ewave = 40              # Energy of packet (eV)
pot_shape = 0           # Shape of potential
                        # 0 = No potential
                        # 1 = Step potential
                        # 2 = Single wall
                        # 3 = Double wall
                        # 4 = Finite well (Packet starts at middle)
pot_height_eV = 20      # eV
barrier_thickness = 10  # Thickness of the barrier(Angstrom = 10^-10 m)
                        # Only for Potential_shape = 2 or 3!
algorithm = 2           # 0 = Forward
                        # 1 = Backward
			# 2 = Crank-Nicolson
