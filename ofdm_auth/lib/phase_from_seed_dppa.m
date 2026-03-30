function phi_deg = phase_from_seed_dppa(seed_u32, phi_choices_deg)
    idx = 1 + mod(double(seed_u32), numel(phi_choices_deg));
    phi_deg = phi_choices_deg(idx);
end