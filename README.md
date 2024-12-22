# QRS Detection

## Overview
This repository contains the implementation of a robust digital QRS-detection algorithm for arrhythmia monitoring, based on the method proposed by [Lindenberg and Kunt](https://www.sciencedirect.com/science/article/pii/0010480983900277). The project was developed as part of the **Biomedical Signal and Image Processing** course at the University of Ljubljana, Faculty of Computer and Information Science, in 2024.

## Project Description
The primary goal of this project is the detection of QRS complexes in ECG signals to accurately extract R peaks. This task is crucial for heartbeat detection and arrhythmia monitoring. The implemented algorithm was tested on:

- [**MIT/BIH Arrhythmia Database**](https://www.physionet.org/content/mitdb/1.0.0/)
- [**Long-Term ST Database**](https://www.physionet.org/content/ltstdb/1.0.0/)

### Key Features
1. **Original Algorithm Implementation**: Based on the robust method by Lindenberg and Kunt.
2. **Weakness Identification**: Analysis of the algorithm's shortcomings.
3. **Proposed Improvements**: Enhancements to the original algorithm to address identified limitations.
4. **Evaluation of Improvements**: Comprehensive testing of both the original and improved algorithms on the aforementioned datasets.

## Repository Structure
- **`Detector.m`**: MATLAB function for testing the QRS detection algorithm on files.
- **`QRSDetect.m`**: MATLAB function implementing the QRS detection algorithm.
- **Report**: A PDF file detailing the algorithm, improvements made, and performance evaluation.
