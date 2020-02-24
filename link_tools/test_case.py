import link_tools

test = {
    'antenna' : link_tools.antenna(),
    'device_name' : "test",
    'tracking' : False,
    'output_power_dBW' : 0,
    'transmission_line_loss' : - 1.0,
    'switch_loss' : - 0.2,
    'pointing_loss' : -0.4,
    'polarization_loss' : - 3.5,
    'radome_loss' : 0,
    'LNA_noise_temp' : 100,
    'transmission_line_temp' : 290,
    'sky_temp' : 275,
    'minimum_rx_SNR' : 15.9,
    'path_loss' : 153,
    'atmospheric_loss' : - 2.1,
    'rain_loss' : 0,
    'ionospheric_loss' : - 0.35,
    'frequency' : 400e6, # Hertz
    'bandwidth' : 40e3,
    'desired_data_rate_bps' : 38400,
    'orbit_name' : 'default',
    'satellite_mass' : 0,
    'apoapsis_height' : 400, # kilometers
    'periapsis_height' : None,
    'body_orbited_name' : 'Earth',
    'path_distance' : 400e3 # meters
}

antenna_test = {
    'peak_gain': 23.5,
    'polarization': 'Linear'
}

device1 = link_tools.rf_device(
                 name = test['device_name'],
                 tracking = test['tracking'],
                 antenna = test['antenna'],
                 output_power_dBW = test['output_power_dBW'],
                 transmission_line_loss = test['transmission_line_loss'],
                 switch_loss = test['switch_loss'],
                 pointing_loss = test['pointing_loss'],
                 polarization_loss = test['polarization_loss'],
                 radome_loss = test['radome_loss'],
                 LNA_noise_temp = test['LNA_noise_temp'],
                 transmission_line_temp = test['transmission_line_temp'],
                 sky_temp = test['sky_temp'],
                 minimum_rx_SNR = test['minimum_rx_SNR']
)

print("Test Device Name: ", device1.name)
print("G over T: ", device1.G_over_T, "dB")
print("Setting new antenna")
device1.set_antenna(link_tools.antenna(
                 peak_gain = antenna_test['peak_gain'],
                 polarization = antenna_test['polarization']))
print("New G over T: ", device1.G_over_T, "dB")
orbit1 = link_tools.satellite_orbit(
                 name = test['orbit_name'],
                 apoapsis_height = test['apoapsis_height'],
                 periapsis_height = test['periapsis_height'],
                 body_orbited_name = test['body_orbited_name'],
                 satellite_mass = test['satellite_mass'])
print("Height of Periapsis: ", orbit1.periapsis_height, "km")
print("Eccentricity of Orbit: ", orbit1.eccentricity)
print("Free Space Path Loss Test: ", link_tools.free_space_path_loss(400e3,400e6), "dB")
test_path = link_tools.path(
                 atmospheric_loss = test['atmospheric_loss'],
                 ionospheric_loss = test['ionospheric_loss'],
                 rain_loss = test['rain_loss'],
                 path_distance = test['path_distance'])
test_message = link_tools.message(desired_data_rate_bps = test['desired_data_rate_bps'],
                frequency = test['frequency'],
                bandwidth = test['bandwidth'])
print("Power density: ", test_path.calculate_power_density(test_message, device1), "W/m^2")
print("EIRP in dBW: ", device1.EIRP_dBW, "dBW")
