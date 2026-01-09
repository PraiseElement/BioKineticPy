# Unit definitions: Factors to convert INPUT to STANDARD (Molar, Seconds)
CONC_TO_MOLAR = {
    "M": 1.0,
    "mM": 1e-3,
    "μM": 1e-6,
    "nM": 1e-9,
    "pM": 1e-12
}

TIME_TO_SECONDS = {
    "sec": 1.0,
    "min": 60.0
}

def get_conversion_factor(conc_unit: str, time_unit: str = "sec"):
    c_factor = CONC_TO_MOLAR.get(conc_unit, 1.0)
    t_factor = TIME_TO_SECONDS.get(time_unit, 1.0)
    return c_factor, t_factor

def convert_velocity_to_standard(v_value, conc_unit, time_unit):
    """v_std (M/s) = v_user * (c_factor / t_factor)"""
    c_factor, t_factor = get_conversion_factor(conc_unit, time_unit)
    return v_value * (c_factor / t_factor)

def convert_param_to_user(value, param_type, conc_unit, time_unit):
    c_factor, t_factor = get_conversion_factor(conc_unit, time_unit)
    
    if param_type == 'concentration':
        return value / c_factor
    elif param_type == 'velocity':
        return value / (c_factor / t_factor)
    return value