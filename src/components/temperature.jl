# -----------------------------------------------------------
# Global Surface Temperature Change
# -----------------------------------------------------------

@defcomp temperature begin

    earth_radius        = Parameter(description="Earth Radius", default=6371000, unit="m") # Earth radius in m
    seconds_per_year    = Parameter(description="Number of seconds in a year", default=60*60*24*365.24219, unit="seconds") # Number of seconds in a year (s yr⁻¹).
    #d   = Parameter(index=[3])      # Thermal response timescales: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (years).
    decay_factor = Parameter(index=[3]) # Thermal response decay factor, calculated as exp(-1/d) where d represents thermal response timescale of jth thermal box.
    q       = Parameter(index=[3])      # Radiative forcing coefficient: [1] Thermal equilibration of deep ocean & [2] Thermal admustment of upper ocean (K W⁻¹m²).
    F       = Parameter(index=[time])   # Total radiative forcing (Wm⁻²).
    ocean_heat_capacity = Parameter(index=[2], description="Heat capacities of Mixed and Deep Layers of the Ocean",
							unit="W/(m^2 y)") # Heat capacities for mixed layer and deep ocean (W m⁻²yr⁻¹).
    Tj_0    = Parameter(index=[3])
    T_0     = Parameter()

    T   = Variable(index=[time])    # Global mean surface temperature anomaly (K).
    Tj  = Variable(index=[time,3])  # Temperature change for three thermal pools (K).
	ntoa_joule = Variable() # TODO: Fill in description with units.
    c_dtemp = Variable(index=[time])   # TODO: Fill in description with units.	
    del_ohc = Variable(index=[time])   # TODO: Fill in description with units.

	# Set up some terms that depend on the sampled ocean heat parameters, but otherwise remain constant in all periods.
    # TODO: Add comments for each calculation.
    function init(p, v, d)
        v.ntoa_joule = 4 * π * p.earth_radius^2 * p.seconds_per_year
	end

    function run_timestep(p, v, d, t)

        if is_first(t)

        	# Set initial condition for three thermal boxes.
        	v.Tj[t,:] = p.Tj_0

        	# Set initial temperature.
        	v.T[t] = p.T_0

        else

        	#Calculate temperature change for the three different thermal response times.
            for j=1:3
                v.Tj[t,j] = v.Tj[t-1,j] * p.decay_factor[j] + p.F[t] * p.q[j] * (1.0 - p.decay_factor[j])
            end
			# TODO: Add calculation comments.
			v.c_dtemp[t] = p.ocean_heat_capacity[1]*(sum(v.Tj[t,:])-sum(v.Tj[t-1,:])) + p.ocean_heat_capacity[2]*(sum(v.Tj[t,:])-sum(v.Tj[t-1,:]))
			# TODO: Add calculation comments.
			v.del_ohc[t]  = v.ntoa_joule * v.c_dtemp[t]	

            #Calculate global mean surface temperature anomaly.
            v.T[t] = sum([v.Tj[t-1,:] v.Tj[t,:]]) / 2
        end
    end
end
