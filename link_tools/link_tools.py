#!/usr/bin/env python
# coding: utf-8

# # Link Tools
# Designed for free space radio (RF) and optical communications
# In progress. Need to be careful about units of inputs


import numpy as np
import math
import re # To be used for accepting more diverse input. Not impelemented at moment

# Constants
gravitational_constant = 6.674 * (10**-11) # N*m^2/kg^2
gravitational_constant_km = gravitational_constant * (10**-9)
boltzmann_constant = -228.6 # dBW/K/Hz

# Utilities for various conversions

def watts_to_dBW(value_in_watts):
    if (value_in_watts == 0):
        return 0
    else:
        return 10.0*math.log(value_in_watts)

def dBW_to_watts(value_in_dBW):
    return 10.0**(value_in_dBW/10.0)

def dBW_to_dBm(value_in_dBW):
    return value_in_dBW + 30.0

def ERP_to_EIRP_dBW(value_in_ERP):
    return value_in_ERP + 2.15

def EIRP_to_ERP_dBW(value_in_EIRP):
    return value_in_EIRP - 2.15

def bps_to_dBHz(value_in_bps):
    return 10.0*math.log10(value_in_bps)

def free_space_path_loss(distance_m, frequency_hz):
    return (20.0*math.log10(distance_m) + 20.0 * math.log10(frequency_hz)
        + 20.0 * math.log10(4*np.pi/3e8))

def frequency_converter(frequency, desired_unit="hertz"):
    # Take the input, convert to desired unit
    # Unimplemented

    prefix_conversion_factors = {
        'yocto': 0.000000000000000000000001,
        'zepto': 0.000000000000000000001,
        'atto': 0.000000000000000001,
        'femto':  0.000000000000001,
        'pico': 0.000000000001,
        'nano': 0.000000001,
        'micro': 0.000001,
        'milli': 0.001,
        'centi': 0.01,
        'deci': 0.1,
        'base': 1.0,
        'deca': 10,
        'hecto': 100,
        'kilo': 1000,
        'mega': 100000,
        'giga': 1000000000,
        'tera': 1000000000000,
        'peta': 1000000000000000,
        'exa' : 1000000000000000000,
        'zetta' : 1000000000000000000000,
        'yotto' : 1000000000000000000000000
    }
    pass

def angle_converter(angle, desired_unit="radians"):
    # Take the input, convert to desired unit
    # Unimplemented
    pass
def distance_converter(distance, desired_unit="meters"):
    # Take the input, convert to desired unit
    # Unimplemented
    pass

bodies_to_orbit = {
    "Earth" : {
        "units" : "kg,km",
        "mass" : 5.9722 * (10**24),
        "radius" : 6378.17
    },
    "Moon" : {
        "units": "kg,km",
        "mass" : 7.342 * (10**22),
        "radius" : 1737.4
    },
    "Mars" : {
        "units": "kg,km",
        "mass" : 6.4171 * (10**23),
        "radius" : 3389.5
    },
    "Venus" : {
        "units": "kg,km",
        "mass" : 4.8675 * (10**23),
        "radius" : 6051.8
    }
}

class antenna():
    """
    Holds information for antennas.

    Inputs:
    * The peak gain of the antenna. Gain is how focused the energy is
        (e.g. fire-hose vs. water balloon)
    * The polarization of the antenna. At the moment the effect of polarization alignment is not
        factored in

    Possible enhancements:
    * Including a radiation pattern parameter

    """

    def __init__(self,
                 peak_gain = None,
                 polarization = None):
        self.antenna_gain = peak_gain if peak_gain is not None else 0 # dB
        self.polarization = polarization if polarization is not None else 'Omni'


class comm_device():
    """
    Base class for communication devices, both optical and RF

    Inputs:
    * Name - Any string to label the device
    * Output Power in dBW - For the calculations, must be given in dBW

    Possible enhancements:
    * Allow for any output power unit

    """
    def __init__(self,
                name = None,
                output_power_dBW = None):
        self.name = name if name is not None else "undefined"
        self.output_power_dBW = output_power_dBW if output_power_dBW is not None else 0


class optical_device(comm_device):
    """
    Optical/Laser devices are generally defined by different terms than their RF
        counterparts. This relies on base characteristics of a comm device.

    Inputs:
    * Name. Passed to comm device super class.
    * Output Power. Passed to comm device super class.
    * Wavelength. Assummed given in nanometers. For reference, visible light is ~400-700nm
    * Beam Divergence. How much does the beam spread. Could use 1/e or 1/e^2
    * Waist diameter. Aka the exit port diameter. Basically how big is the hole the light
        comes out of

    Internal methods:
    * Calculate power density. Done at a certain distance. As the beam travels, it'll spread out
        by the beam divergence. Farther from the exit port, less power density. This is key for
        laser safety calculations. Assumeds free-space and ignores cloud, etc. attenuations

    * Calculate spot size. Done at a certain distance. Utilized in power density calculations.
        Beam will spread over distance, this calculates the size of that spread at a certain
        distance.

    """

    def __init__(self,
                name = None,
                output_power_dBW = None,
                wavelength = None,
                beam_divergence = None, # In radians
                waist_diameter = None
                ):
        super().__init__(name, output_power_dBW)
        self.wavelength = wavelength if wavelength is not None else 1550e-9 # nanometers
        self.beam_divergence = beam_divergence if beam_divergence is not None else 0
        self.waist_diameter = waist_diameter if waist_diameter is not None else 0

    def calculate_power_density(self, distance_in_km):
        spot_size = calculate_spot_size(distance_in_km)
        power_density = dBW_to_watts(self.output_power_dBW)/spot_size
        return power_density

    def calculate_spot_size(self, distance_in_km):
        # Spot diameter uses Diameter_at_target = diameter_at_exit + 2*
        # (distance_to_target*tangent(beam_divergence/2))
        spot_diameter = self.waist_diameter + 2*(distance_in_km*1000 * math.tan(self.beam_divergence/2))
        spot_size = np.pi*((spot_diameter/2)**2)
        return spot_size


class rf_device(comm_device):
    """
    RF devices are generally defined by different terms than their laser/optical
        counterparts. This relies on base characteristics of a comm device.

    Inputs:
    * Name. Passed to comm device super class.
    * Output Power. Passed to comm device super class.
    * Antenna. Pass in an instance of the antenna class.
    * Tracking. Does the device try to keep pointing towards a receiver during communication
    * Losses:
        Transmission Line Loss. How much noise is there in the transmission line (in dB)
        Switch Loss.
        Pointing Loss. Any inherent/expected inaccuracy in the pointing of the antenna
        Polarization Loss. Can be calculated based on receiver+transmitter polarization.
            Assumed it is given in dB at this point.
        Radome Loss.
    * Temperatures: (Higher temp indicates higher noise floor)
        LNA Noise Temp
        Transmission Line Temp
        Sky Temp
    * Minimum SNR. SNR = signal-to-noise. Basically how much stronger does the sent signal
        need to be than the noise floor. This is generally set by signal processing scheme
        and rule-of-thumb buffer.

    Internal methods:
    * Calculate characertistics. There are a few characteristics that are caclulated from inputs
        This method captures all of those. Performed at initialization. Can be called again
        if inputs change.

    * Set antenna. Can be used to swap out antennas

    """
    def __init__(self,
                 name = None,
                 output_power_dBW = None,
                 antenna = None,
                 tracking = None,
                 transmission_line_loss = None,
                 switch_loss = None,
                 pointing_loss = None,
                 polarization_loss = None,
                 radome_loss = None,
                 LNA_noise_temp = None,
                 transmission_line_temp = None,
                 sky_temp = None,
                 minimum_rx_SNR = None
                ):

        super().__init__(name, output_power_dBW)

        # Characteristics
        self.antenna = antenna
        self.tracking = tracking if tracking is not None else False # Does the device track its target via antenna or device pointing

        self.transmission_line_loss = transmission_line_loss if transmission_line_loss is not None else 0
        self.switch_loss = switch_loss if switch_loss is not None else 0

        self.pointing_loss = pointing_loss if pointing_loss is not None else 0
        self.polarization_loss = polarization_loss if polarization_loss is not None else 0
        self.radome_loss = radome_loss if radome_loss is not None else 0
        self.LNA_noise_temp = LNA_noise_temp if LNA_noise_temp is not None else 100
        self.transmission_line_temp = transmission_line_temp if transmission_line_temp is not None else 290
        self.sky_temp = sky_temp if sky_temp is not None else 275
        self.minimum_rx_SNR = minimum_rx_SNR if minimum_rx_SNR is not None else 1

        self.calculate_characteristics()

    def calculate_characteristics(self):
        self.transmission_line_coeff = 10**(self.transmission_line_loss/10)
        self.effective_noise_temp = ((self.transmission_line_coeff * self.sky_temp) +
                                     (1 - self.transmission_line_coeff) * self.transmission_line_temp +
                                     self.LNA_noise_temp)
        if (self.antenna == None):
            self.EIRP_dBW = None
            self.G_over_T = None
            print("Warning: No antenna attached. Antenna dependent characteristics left undefined.")
        else:
            self.G_over_T = (self.radome_loss + self.antenna.antenna_gain + self.transmission_line_loss -
                            10*math.log10(self.effective_noise_temp))
            self.EIRP_dBW = (self.output_power_dBW + self.transmission_line_loss +
                             self.switch_loss + self.antenna.antenna_gain)

    def set_antenna(self, antenna):
        self.antenna = antenna
        self.calculate_characteristics()


class message():
    """
    This refers to the proposed content of the signal. For the moment, this is soley for RF

    Inputs:
    * Desired Data Rate BPS. In bits-per-second. The higher the BPS, the more
        sensitive the transmission
    * Frequency. This is the electromagnetic frequency of the message, e.g. 400 MHz
    * Bandwidth. Essentially, how much frequency can the signal be spread out over.

    Internal methods:
    * Calculate characertistics. This is a general method name used to caclulate certain things
        based off inputs.
    """

    def __init__(self,
                desired_data_rate_bps = None,
                frequency = None,
                bandwidth = None):
        self.desired_data_rate_bps = desired_data_rate_bps if desired_data_rate_bps is not None else 0
        self.frequency = frequency if frequency is not None else 0
        self.bandwidth = bandwidth if bandwidth is not None else 0

        self.calculate_characteristics()

    def calculate_characteristics(self):
        self.desired_data_rate_dBHz = bps_to_dBHz(self.desired_data_rate_bps)


class path():
    """
    This refers to the proposed path of the signal. For example, sending the signal through
        a cloudy sky vs. free space will lead to a lot more loss

    Inputs:
    * Atmospheric loss. More clouds or turbulence leads to more loss
    * Ionospheric loss.
    * Rain loss. If there is rain, this will lead to losses. Generally assume a certain level
    * Path distance. How far is the signal travelling

    Internal methods:
    * Calculate budget. This calculates the link budget based on a given device transmitting,
        receiving device, message, and path.
    * Calculate power density. Note this does the equivalent to the optical power density but
        bases spread on idealized omni-directional spread as opposed to on beam divergence

    Possible enhancements:
    * Be able to input a date, time, place or some other parameter to then pull associated losses
    * Link this closer to satellite orbit so that path distance varies over course of a pass
    * Make more applicable to both optical and RF
    """

    def __init__(self,
                 atmospheric_loss = None,
                 ionospheric_loss = None,
                 rain_loss = None,
                 path_distance = None):
        self.atmospheric_loss = atmospheric_loss if atmospheric_loss is not None else 0
        self.ionospheric_loss = ionospheric_loss if ionospheric_loss is not None else 0
        self.rain_loss = rain_loss if rain_loss is not None else 0
        self.path_distance = path_distance if path_distance is not None else 0

    def calculate_budget(self, transmit_device, receive_device, message):
        self.path_loss = free_space_path_loss(self.path_distance, message.frequency)
        self.isotropic_signal_at_rx = (transmit_device.EIRP_dBW +
                                       transmit_device.pointing_loss +
                                       transmit_device.polarization_loss +
                                       self.path_loss +
                                       self.atmospheric_loss +
                                       self.ionospheric_loss +
                                       self.rain_loss)
        self.rx_signal_to_noise_power_density = (self.isotropic_signal_at_rx +
                                                receive_device.pointing_loss +
                                                -boltzmann_constant +
                                                receive_device.G_over_T)
        self.power_at_LNA_receive = (isotropic_signal_at_rx +
                                receive_device.pointing_loss +
                                receive_device.radome_loss +
                                receive_device.antenna_gain +
                                receive_device.transmission_line_loss)
        self.receive_noise_power = (boltzmann_constant +
                                    10*math.log10(receive_device.effective_noise_temp) +
                                   10*math.log10(message.bandwidth))
        self.budget_SNR = self.rx_signal_to_noise_power_density - message.desired_data_rate_dBHz
        self.link_margin = self.receive_noise_power - receive_device.minimum_rx_SNR
        return self.link_margin

    def calculate_power_density(self, message, transmit_device):
        self.path_loss = free_space_path_loss(self.path_distance, message.frequency)
        self.power_density = ((dBW_to_watts(transmit_device.EIRP_dBW)) /
                             (4.0 * np.pi * (self.path_distance**2)))
        return self.power_density


class satellite_orbit():
    """
    This refers to orbit of a satellite. Able to be used for any inputted
        body parameters (e.g. Earth, Moon, Mars)

    Inputs:
    * Name. Used to provide some label for the orbit (e.g. ISS Orbit)
    * Apoapsis Height. What's the farthest the satellite is from what it's orbiting.
        Note: The Apo- prefix refers to farthest point. The suffix generally refers to the body
        oribited. E.g. apogee is farthest point in an Earth orbit. -apsis suffix is general
    * Periapsis Height. What's the closest the satellite is from what it's orbiting.

    Internal methods:
    * Calculate characteristics. Used to calculate a few parameters based on inputs. Called
        during initialization. Can be called again with any updates.

    Possible enhancements:
    * Make general for a wider variety of orbits
    """
    def __init__(self,
                 name = None,
                 apoapsis_height = None,
                 periapsis_height = None,
                 body_orbited_name = None,
                 satellite_mass = None):

        self.name = name if name is not None else "undefined"
        self.body_orbited_name = body_orbited_name if body_orbited_name is not None else "Earth"
        self.body_orbited = bodies_to_orbit[self.body_orbited_name]
        self.apoapsis_height = apoapsis_height if apoapsis_height is not None else 0
        self.satellite_mass = satellite_mass if satellite_mass is not None else 0
        self.periapsis_height = periapsis_height if periapsis_height is not None else self.apoapsis_height

        self.calculate_characteristics()

    def calculate_characteristics(self):
        self.semi_major_axis = self.apoapsis_height + self.body_orbited['radius']
        self.semi_minor_axis = self.periapsis_height + self.body_orbited['radius']
        self.mean_altitude = (self.apoapsis_height+self.periapsis_height)/2
        self.orbital_period = (2*np.pi) * np.sqrt(self.semi_major_axis**3/
                                                  (gravitational_constant_km*(self.body_orbited["mass"]+self.satellite_mass)))
        self.eccentricity = (self.semi_major_axis-self.semi_minor_axis)/(self.semi_major_axis+self.semi_minor_axis)
