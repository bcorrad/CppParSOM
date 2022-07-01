#ifndef UNTITLED_GLOBALS_H
#define UNTITLED_GLOBALS_H

constexpr int MAX_THREADS = 4;
constexpr int FLOAT_MIN = 0;
constexpr int FLOAT_MAX = 1;
constexpr int WEIGHT_SIZE = 100; // size of weight vector of each node [10, 100, 1000].
constexpr int EPOCHS = 10;   //the number of epochs desired for the training [10, 100, 1000].
constexpr double START_LR = 0.001;    //the value of the learning rate at the start of training.
constexpr int GRID_ROWS = 100;
constexpr int GRID_COLS = 100;
constexpr const int N_INPUTS = 1000; // number of weights each node must contain. [100, 1000, 10000]

#endif //UNTITLED_GLOBALS_H
