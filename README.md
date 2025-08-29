# TI-Stim-Waveform-Generation
MATLAB code for generating waveforms and controlling stimulation/imaging hardware for **temporal interference (TI) stimulation** experiments. This repository supports experiments such as those described in [Song et al., *Nature Neuroscience*, 2023](https://www.nature.com/articles/s41593-023-01456-8).

## Overview
The primary purposes of this code are:
- Generate temporal interference stimulation waveforms (pulse trains).  
- Control connected hardware (e.g., NI-DAQ, imaging laser, camera trigger).  
- Log experiment status and optionally visualize waveform summaries.  
The repository is designed to run with minimal dependencies:  
only MATLAB R2018a – R2025a is required.

## Features
- Configurable pulse train generation with ramped onset/offset  
- Flexible amplitude and frequency settings  
- Control of NI-DAQ outputs for:
  - Stimulation signals  
  - Imaging laser trigger  
  - External device triggers (e.g., camera, electrophysiology hardware)  
- Built-in logging and plotting utilities  
- Simple parameter editing via a single script  

## Requirements
- MATLAB R2018a – R2025a  
- NI-DAQ device (tested with NI-DAQmx devices)  
No additional toolboxes or packages are required.

## Quick Demo
The easiest way to try the code is to run the included script:
```matlab
% Run a demo TI stimulation experiment (generates waveforms + plots)
RUNME_pulse_train_ti_stimulation_experiment
```
This script will:
- Generate a pulse train waveform pair
- Plot the individual signals and their summed interference
- (If hardware is connected) send the signals to the NI-DAQ device
**Note:** To configure parameters (pulse duration, amplitudes, frequencies, etc.),
open `source/RUNME_pulse_train_ti_stimulation_experiment.m` in MATLAB and edit the
`inputSpecification` section directly before running.

### Example: Updating Parameters
To increase the pulse frequency from 1 Hz to 5 Hz,
find the inputSpecification section and update:

```matlab
% The frequency of pulses in the stimulation (Hz)
'PulseFreq', 1   % original
```
to
```matlab
'PulseFreq', 5   % updated
```
Re-run the script to generate new waveforms and see the changes reflected in the plots.

## Main Function: `createTIPulseTrainWaveforms`
The core waveform generator function is:
```matlab
[waveforms, t, params] = createTIPulseTrainWaveforms(paramStruct)
```
where `paramStruct` is a MATLAB struct of experiment parameters.
The most important parameters are:
| Parameter              | Type             | Units | Description                                                                    |
| ---------------------- | ---------------- | ----- | ------------------------------------------------------------------------------ |
| `PulseDuration`        | scalar           | s     | Duration of each stimulation pulse. Must be `< 1/PulseFreq`.                   |
| `RampDuration`         | scalar           | s     | Duration of ramp-up and ramp-down at start and end of pulses.                  |
| `PulseFreq`            | scalar           | Hz    | Frequency of pulses within a train (not the TI beat frequency).                |
| `PulsesPerTrain`       | scalar           | –     | Number of pulses per train.                                                    |
| `CarrierFreq`          | scalar           | Hz    | Carrier frequency of the TI waveforms.                                         |
| `InterferenceBeatFreq` | scalar           | Hz    | Frequency of the interference envelope.                                        |
| `NumTrains`            | scalar           | –     | Number of pulse trains to run in one experiment.                               |
| `InterTrainInterval`   | scalar           | s     | Interval between pulse trains.                                                 |
| `A1`, `A2`             | scalar or vector | V     | Amplitudes of the two signals. Can be vectors (one per train).                 |
| `SamplingRate`         | scalar           | Hz    | Sampling rate of the generated signals.                                        |
| `Flip`                 | boolean          | –     | If true, flips one waveform by 180° (phase correction).                        |
| `WaitTime`             | scalar           | s     | Pre/post-stimulation wait time.                                                |
| `ModulationSymetry`    | string           | –     | `'signal1'`, `'signal2'`, or `'symetric'`. Controls how modulation is applied. |
| `Control`              | boolean          | –     | Generate an unmodulated control waveform.                                      |
| `Plot`                 | boolean          | –     | Plot waveform summaries if true.                                               |
| `Debug`                | boolean          | –     | Enable debug mode with extra logging.                                          |

The function returns:
- waveforms: generated stimulation waveforms (matrix, samples × channels)
- t: time vector for the generated signals
- params: parameter struct actually used in generation

## Output
Running an experiment produces:
- Command line logs of experiment status
- Optional summary plots of generated waveforms
- Live stimulation signals delivered through NI-DAQ (if connected)

## Contributing
Feature requests, bug reports, and hardware test results are welcome.
Please open an issue or pull request on GitHub.

## License
This project is provided under the MIT License.

## Attributions
This code was developed and written by Corban Swain. `README.md` file written with the assistance of ChatGPT (GPT-5).
