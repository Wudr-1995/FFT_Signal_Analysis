# FFT_Signal_Analysis
A signal processing method using FFT.

## Compile Need
CERN ROOT & JUNO offline

##Compile Command
```bash
make
```

And clean compile

```bash
make clean
```

## Command
1. "./getParas".
2. Step (1 or 2).
3. Input type (1 or 2).
4. Input file name.
5. Lower bound.
6. Higher bound.
7. Scale.
8. The number of sampling points.
9. The number of waveforms.
10. Amplitude threshold.
11. Charge threshold.
12. Baseline input file name (just for step 1).

### Input file

Waveforms.txt file (this file content signal and baseline waveforms)

Baseline.txt file (this file just content baseline waveforms)

Because of the data requirment system we use, there are always two column in the input txt, the first column is the data we need, the other is useless. Just like this.

-0.001875	-0.001875

-0.002031	-0.002031

-0.002188	-0.002188

-0.002031	-0.002031

-0.001875	-0.001875

-0.001875	-0.001875

-0.001719	-0.001719

-0.002031	-0.002031

-0.002500	-0.002500

-0.002188	-0.002188

-0.001719	-0.001719

-0.001719	-0.001719

-0.002500	-0.002500

The program also can accept the input file who has only one column. It's necessary to set the first argument input in the command line (1: one line type, 2: two line type). If type 1 is chosen, the program also can due with input file who has one waveform per line.

For convenience, all waveform data were put into one txt file. (This step can be achieved by "type i\*.txt >> ../input.txt" command in windows cmd")

### An example
```bash
./getParas 1 2 1m_1305 400 800 0.2 1400 10000 2.5 0.2 baseline
./getParas 2 2 1m_1305 400 800 0.2 1400 10000 1 0.2
```
The first step is to get a Wiener filter for these signals.
The second step is to use the filter we get in last step to due with signals and get parameters of signals.
