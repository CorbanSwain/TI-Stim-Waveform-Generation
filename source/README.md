# Temporal Interference (TI) Stimulation Waveform Generation and Experiment Control

This repository contains a MATLAB script for configuring stimulation, generating waveforms, and running Temporal Interference (TI) stimulation experiments.

## Getting Started

To use this script, follow these steps:

1. Clone or download the repository from GitHub.

2. Ensure you have MATLAB installed on your computer.

3. Open MATLAB and navigate to the `source` directory where the script is located.

4. Run the script by executing the entire code or running specific sections as needed.

## Parameters for Running

The script accepts various parameters to customize the stimulation and experiment settings. Here is a brief description of each parameter:

- `Duration`: The duration of the stimulation in seconds.
- `CarrierFreq`: The carrier frequency of the signals in Hz.
- `BreakCarrierFreq`: The carrier frequency during the break of the signals in Hz (or 'same' to use the same as `CarrierFreq`).
- `PulseFreq`: The frequency of the pulses in Hz.
- `DutyCycle`: The duty cycle of the pulses. A value of 1 means continuous pulses, while 0 means no pulses are present.
- `A1`: The amplitude of the first signal, S1, in volts.
- `A2`: The amplitude of the second signal, S2, in volts.
- `SamplingRate`: The sampling rate in Hz used for generating the signals.
- `RampTime`: The duration of the ramp-up and ramp-down period of the signals in seconds.
- `Flip`: Whether or not to flip one of the signals (shift in phase 180 deg) - true or false.
- `FMTime`: The modulation time during which the waveforms will smoothly transition from an initial frequency and relative phase to a new frequency and relative phase in seconds (or 'auto').
- `WaitTime`: The duration of time before and after the stimulation starts (pre-stimulation time) in seconds.
- `ModulationSymmetry`: Controls which signals are modulated to produce the stimulation pulses. Valid values are 'signal1', 'signal2', or 'symmetric'.
- `Control`: A flag to keep the two waveforms out of phase for the entire stimulation - true or false.
- `Plot`: A flag to plot a summary of the output signals - true or false.
- `FMArgs`: A cell array of Name, Value arguments sent to `generateFMPhaseFcn`.
- `Debug`: A flag to run the code in debug mode - true or false.

## Required Hardware

To fully execute the code, you will need the following hardware:

- National Instruments Data Acquisition (DAQ) device (e.g., Dev2).
- Imaging laser and its control interface.

## Resources for Debugging and Submitting Issues

If you encounter any issues or need help with debugging, you can refer to the following resources:

1. GitHub Issues: Visit the [repository's Issues section](/issues) to check for existing issues or submit a new one.

2. Documentation: You can find additional information and explanations in the source code comments and MATLAB documentation.

3. Contact: If you need further assistance or have specific questions, you can contact the repository owner (Corban Swain).

## Acknowledgments

This README.md file was created with the help of Open AI's ChatGPT.
