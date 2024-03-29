#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

/*
Algorytm jest ten sam, więc bardziej szczegółowe wyjaśnienia są w pliku seq.cpp
(komentarze) Różnią się częsci z symulacja
*/

#define C sqrt(2)
#define SIZE 9
#define NUM_OF_MOVEMENTS_IN_SIMULATION 10
#define MAX_DEPTH 100
#define MOVEMENTS 40
#define MAX_NUMBER_OF_THREADS 1024

enum State { EMPTY, BLACK, WHITE };

std::vector<std::pair<int, int>> NEIGHBOURS[SIZE][SIZE];
State previousPositionForBlack[SIZE][SIZE];
State previousPositionForWhite[SIZE][SIZE];
std::chrono::duration<double> total_time_selection;
std::chrono::duration<double> total_time_expansion;
std::chrono::duration<double> total_time_simulation;
std::chrono::duration<double> total_time_backpropagation;
std::chrono::duration<double> total_time_CPU_GPU_copying;
long selection_moves = 0;
long expansion_moves = 0;
long simulation_moves = 0;
long backpropagation_moves = 0;
int total_taken_black = 0;
int total_taken_white = 0;

typedef struct Node {
  State board[SIZE][SIZE];
  std::vector<Node *> children = std::vector<Node *>();
  Node *parent = NULL;
  int taken_black_stones = 0;
  int taken_white_stones = 0;
  unsigned int number_of_simulations = 0;
  double black_score = 0.0;
  ~Node() {
    for (Node *child : children) {
      delete child;
    }
  }

  std::vector<Node *> getMaxProfitForWhite() {
    std::vector<Node *> profitChildren;
    int profit = taken_black_stones - taken_white_stones;
    int tmp = 0;
    for (Node *child : children) {
      tmp = child->taken_black_stones - child->taken_white_stones;
      if (tmp > profit) {
        profit = tmp;
        profitChildren.clear();
        profitChildren.push_back(child);
      } else if (tmp == profit) {
        profitChildren.push_back(child);
      }
    }
    return profitChildren;
  }

  std::vector<Node *> getMaxProfitForBlack() {
    std::vector<Node *> profitChildren;
    int profit = taken_white_stones - taken_black_stones;
    int tmp = 0;
    for (Node *child : children) {
      tmp = child->taken_white_stones - child->taken_black_stones;
      if (tmp > profit) {
        profit = tmp;
        profitChildren.clear();
        profitChildren.push_back(child);
      } else if (tmp == profit) {
        profitChildren.push_back(child);
      }
    }
    return profitChildren;
  }
} Node;

double calculateUct(Node *n, State state) {
  if (n->number_of_simulations == 0) {
    return std::numeric_limits<double>::infinity();
  }
  double w_i = 0.0;
  if (state == BLACK) {
    w_i = n->black_score;
  } else {
    w_i = (double)n->number_of_simulations - n->black_score;
  }
  return (w_i / n->number_of_simulations) +
         C * sqrt(log(n->parent->number_of_simulations) /
                  n->number_of_simulations);
}

void copyBoard(State source[SIZE][SIZE], State destination[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      destination[i][j] = source[i][j];
    }
  }
}

void printBoard(State board[SIZE][SIZE]) {
  std::cout << "    ";
  for (int i = 0; i < SIZE; ++i) {
    std::cout << i << ' ';
  }
  std::cout << '\n';
  for (int i = 0; i < SIZE; ++i) {
    std::cout << i << ' ' << ' ' << " ";
    for (int j = 0; j < SIZE; ++j) {
      if (board[i][j] == EMPTY) {
        std::cout << ". ";
      } else if (board[i][j] == BLACK) {
        std::cout << "X ";
      } else {
        std::cout << "O ";
      }
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

bool isKo(State prev[SIZE][SIZE], State actual[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (prev[i][j] != actual[i][j])
        return false;
    }
  }
  return true;
}

void createNeighbours() {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      std::vector<std::pair<int, int>> neighbours;
      if (i - 1 >= 0) {
        neighbours.push_back(std::make_pair(i - 1, j));
      }
      if (i + 1 < SIZE) {
        neighbours.push_back(std::make_pair(i + 1, j));
      }
      if (j - 1 >= 0) {
        neighbours.push_back(std::make_pair(i, j - 1));
      }
      if (j + 1 < SIZE) {
        neighbours.push_back(std::make_pair(i, j + 1));
      }
      NEIGHBOURS[i][j] = neighbours;
    }
  }
}

std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
findReached(State board[SIZE][SIZE], int i, int j) {
  State color = board[i][j];
  std::vector<std::pair<int, int>> chain, reached;
  std::vector<std::pair<int, int>> frontier = {std::make_pair(i, j)};
  while (!frontier.empty()) {
    auto current_fc = frontier.back();
    frontier.pop_back();
    if (std::find(chain.begin(), chain.end(), current_fc) == chain.end()) {
      chain.push_back(current_fc);
    }
    for (auto fn : NEIGHBOURS[current_fc.first][current_fc.second]) {
      if (board[fn.first][fn.second] == color &&
          std::find(chain.begin(), chain.end(), fn) == chain.end()) {
        frontier.push_back(fn);
      } else if (board[fn.first][fn.second] != color &&
                 std::find(reached.begin(), reached.end(), fn) ==
                     reached.end()) {
        reached.push_back(fn);
      }
    }
  }
  return std::make_pair(chain, reached);
}

State changeState(State state) {
  if (state == BLACK) {
    return WHITE;
  } else {
    return BLACK;
  }
}

std::pair<bool, std::vector<std::pair<int, int>>>
couldPlaceStone(State board[SIZE][SIZE], int row, int col, State state) {
  std::vector<std::pair<int, int>> taken_stones =
      std::vector<std::pair<int, int>>();
  if (board[row][col] != EMPTY)
    return std::make_pair(false, taken_stones);
  board[row][col] = state;
  State alternativeState = changeState(state);

  std::vector<std::pair<int, int>> stones_reached_by_mine =
      findReached(board, row, col).second;
  // tu sprawdzamy czy cos zabieramy u przeciwnika
  for (int i = 0; i < stones_reached_by_mine.size(); ++i) {
    if (board[stones_reached_by_mine[i].first]
             [stones_reached_by_mine[i].second] !=
        EMPTY) { // dla kazego kamienia z sąsiadujących
      auto tmp_chain_reached =
          findReached(board, stones_reached_by_mine[i].first,
                      stones_reached_by_mine[i].second);
      auto potential_captured = tmp_chain_reached.first;
      auto reached_by_opponent = tmp_chain_reached.second;
      int j = 0;
      for (j = 0; j < reached_by_opponent.size();
           ++j) { // dla kazdego wezla co jest obok kamienia przeciwnika
        if (board[reached_by_opponent[j].first]
                 [reached_by_opponent[j].second] != state) {
          break;
        }
      }
      if (j == reached_by_opponent.size()) {
        for (int ind = 0; ind < potential_captured.size(); ++ind) {
          taken_stones.push_back(std::make_pair(
              potential_captured[ind].first, potential_captured[ind].second));
        }
      }
    }
  }
  // end
  // jak nic nie zabieramy to czy nie ma samobojstwa
  if (taken_stones.size() == 0) {
    int j = 0;
    for (j = 0; j < stones_reached_by_mine.size(); ++j) {
      if (board[stones_reached_by_mine[j].first]
               [stones_reached_by_mine[j].second] != alternativeState) {
        break;
      }
    }
    if (j == stones_reached_by_mine
                 .size()) { // only opponent stones surround me, suicide
      board[row][col] = EMPTY;
      return std::make_pair(false, taken_stones);
    }
  }
  // KO check
  bool is_ko;
  for (int i = 0; i < taken_stones.size(); ++i) {
    board[taken_stones[i].first][taken_stones[i].second] = EMPTY;
  }
  if (state == WHITE) {
    is_ko = isKo(previousPositionForBlack, board);
  } else {
    is_ko = isKo(previousPositionForWhite, board);
  }
  if (is_ko) {
    return std::make_pair(false, taken_stones);
  }

  for (int i = 0; i < taken_stones.size(); ++i) {
    board[taken_stones[i].first][taken_stones[i].second] = alternativeState;
  }
  board[row][col] = EMPTY;
  return std::make_pair(true, taken_stones);
}

void createChildren(Node *n, int i, int j, State state) {
  State tmp_board[SIZE][SIZE];
  copyBoard(n->board, tmp_board);
  auto taken_stones = couldPlaceStone(tmp_board, i, j, state).second;
  Node *childNode = new Node();
  copyBoard(n->board, childNode->board);
  childNode->board[i][j] = state;
  childNode->taken_black_stones = n->taken_black_stones;
  childNode->taken_white_stones = n->taken_white_stones;
  for (int i = 0; i < taken_stones.size(); ++i) {
    childNode->board[taken_stones[i].first][taken_stones[i].second] = EMPTY;
  }
  if (state == BLACK) {
    childNode->taken_white_stones =
        childNode->taken_white_stones + taken_stones.size();
  } else {
    childNode->taken_black_stones =
        childNode->taken_black_stones + taken_stones.size();
  }
  n->children.push_back(childNode);
  childNode->parent = n;
}

std::pair<int, int> computeTerritories(State board[SIZE][SIZE]) {
  bool managed[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      managed[i][j] = false;
    }
  }
  std::vector<std::pair<int, int>> chain;
  std::vector<std::pair<int, int>> reached;
  State color;
  int white_territory = 0;
  int black_territory = 0;
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (board[i][j] == EMPTY && !managed[i][j]) {
        auto res = findReached(board, i, j);
        chain = res.first;
        reached = res.second;
        for (auto p : chain) {
          managed[p.first][p.second] = true;
        }
        // jesli w reached wszystkie kamienie sa jednego koloru to to jest
        // terytorium tego koloru
        color = board[reached[0].first][reached[0].second];
        int k = 0;
        for (k = 0; k < reached.size(); ++k) {
          if (board[reached[k].first][reached[k].second] != color) {
            break;
          }
        }
        if (k == reached.size()) {
          if (color == BLACK) {
            black_territory += chain.size();
          } else {
            white_territory += chain.size();
          }
        }
      }
    }
  }
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (board[i][j] == BLACK) {
        ++black_territory;
      } else if (board[i][j] == WHITE) {
        ++white_territory;
      }
    }
  }
  return std::make_pair(black_territory, white_territory);
}

void expand(Node *n, State state) {
  State tmp_board[SIZE][SIZE];
  copyBoard(n->board, tmp_board);
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (couldPlaceStone(tmp_board, i, j, state).first) {
        createChildren(
            n, i, j,
            state); // i j - tu ustawimy kamien i takie dziecko dodamy do n
      }
    }
  }
}

void flattenCube(Node *n, State *flattenedCube) {
  int index = 0;
  for (int k = 0; k < n->children.size();
       ++k) { // tworzymy kostkę żeby potem przesłać plansze wszystkich dzieci
              // na GPU
    for (int i = 0; i < SIZE; ++i) {
      for (int j = 0; j < SIZE; ++j) {
        flattenedCube[index] = n->children[k]->board[i][j];
        ++index;
      }
    }
  }
}

// analogia findReached, tylko zamiast wektora par używam 2-wymiarowych tablic
// na wejsciu chain i reached są wypełnione -1 (dla łatwiejszego pilnowania
// rozmiaru)
__device__ void d_findReached(State board[SIZE][SIZE], int i, int j,
                              int chain[SIZE * SIZE][2],
                              int reached[SIZE * SIZE][2]) {
  State color = board[i][j];
  int frontier[SIZE * SIZE][2] = {{i, j}};
  int frontierSize = 1;
  int current_fc[2];
  bool in_chain[SIZE][SIZE];
  bool in_reached[SIZE][SIZE];
  int chain_size = 0;
  int reached_size = 0;
  for (int i = 0; i < SIZE * SIZE; ++i) {
    in_chain[i][0] = false;
    in_chain[i][1] = false;
    in_reached[i][0] = false;
    in_reached[i][1] = false;
  }
  while (frontierSize != 0) {
    current_fc[0] = frontier[frontierSize - 1][0];
    current_fc[1] = frontier[frontierSize - 1][1];
    --frontierSize;
    if (!in_chain[current_fc[0]][current_fc[1]]) {
      chain[chain_size][0] = current_fc[0];
      chain[chain_size][1] = current_fc[1];
      in_chain[current_fc[0]][current_fc[1]] = true;
      ++chain_size;
    }
    int i = current_fc[0];
    int j = current_fc[1];

    if (i - 1 >= 0) {
      if (board[i - 1][j] == color && !in_chain[i - 1][j]) {
        frontier[frontierSize][0] = i - 1;
        frontier[frontierSize][1] = j;
        ++frontierSize;
      } else if (board[i - 1][j] != color && in_reached[i - 1][j]) {
        reached[reached_size][0] = i - 1;
        reached[reached_size][1] = j;
        in_reached[i - 1][j] = true;
        ++reached_size;
      }
    }
    if (i + 1 < SIZE) {
      if (board[i + 1][j] == color && !in_chain[i + 1][j]) {
        frontier[frontierSize][0] = i + 1;
        frontier[frontierSize][1] = j;
        ++frontierSize;
      } else if (board[i + 1][j] != color && in_reached[i + 1][j]) {
        reached[reached_size][0] = i + 1;
        reached[reached_size][1] = j;
        in_reached[i + 1][j] = true;
        ++reached_size;
      }
    }
    if (j - 1 >= 0) {
      if (board[i][j - 1] == color && !in_chain[i][j - 1]) {
        frontier[frontierSize][0] = i;
        frontier[frontierSize][1] = j - 1;
        ++frontierSize;
      } else if (board[i][j - 1] != color && in_reached[i][j - 1]) {
        reached[reached_size][0] = i;
        reached[reached_size][1] = j - 1;
        in_reached[i][j - 1] = true;
        ++reached_size;
      }
    }
    if (j + 1 < SIZE) {
      if (board[i][j + 1] == color && !in_chain[i][j + 1]) {
        frontier[frontierSize][0] = i;
        frontier[frontierSize][1] = j + 1;
        ++frontierSize;
      } else if (board[i][j + 1] != color && in_reached[i][j + 1]) {
        reached[reached_size][0] = i;
        reached[reached_size][1] = j + 1;
        in_reached[i][j + 1] = true;
        ++reached_size;
      }
    }
  }
}

// taken_stones jest wypełnione -1
__device__ bool d_couldPlaceStone(State board[SIZE][SIZE], int row, int col,
                                  State state,
                                  int taken_stones[SIZE * SIZE][2]) {
  int taken_stones_size = 0;
  if (board[row][col] != EMPTY) {
    return false;
  }
  board[row][col] = state;
  State alternativeState = BLACK;
  if (state == BLACK) {
    alternativeState = WHITE;
  }

  int chain[SIZE * SIZE][2];
  int reached[SIZE * SIZE][2];
  for (int i = 0; i < SIZE * SIZE; ++i) {
    chain[i][0] = -1;
    chain[i][1] = -1;
    reached[i][0] = -1;
    reached[i][1] = -1;
  }

  d_findReached(board, row, col, chain, reached);
  int r_ind = 0;
  while (reached[r_ind][0] != -1) {
    int potential_captured[SIZE * SIZE][2];
    int reached_by_opponent[SIZE * SIZE][2];
    for (int i = 0; i < SIZE * SIZE; ++i) {
      potential_captured[i][0] = -1;
      potential_captured[i][1] = -1;
      reached_by_opponent[i][0] = -1;
      reached_by_opponent[i][1] = -1;
    }

    if (board[reached[r_ind][0]][reached[r_ind][1]] != EMPTY) {
      d_findReached(board, reached[r_ind][0], reached[r_ind][1],
                    potential_captured, reached_by_opponent);
      int j = 0;
      while (reached_by_opponent[j][0] != -1) {
        if (board[reached_by_opponent[j][0]][reached_by_opponent[j][1]] !=
            state) {
          break;
        }
        ++j;
      }
      if (reached_by_opponent[j][0] == -1) {
        int ind = 0;
        while (potential_captured[ind][0] != -1) {
          taken_stones[taken_stones_size][0] = potential_captured[ind][0];
          taken_stones[taken_stones_size][1] = potential_captured[ind][1];
          ++taken_stones_size;
          ++ind;
        }
      }
    }
    ++r_ind;
  }

  if (taken_stones[0][0] == -1) {
    int j = 0;
    while (reached[j][0] != -1) {
      if (board[reached[j][0]][reached[j][1]] != alternativeState) {
        break;
      }
      ++j;
    }
    if (reached[j][0] == -1) {
      board[row][col] = EMPTY;
      return false;
    }
  }
  board[row][col] = EMPTY;
  return true;
}

__device__ void d_computeTerritories(State board[SIZE][SIZE], int results[2]) {
  bool managed[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      managed[i][j] = false;
    }
  }
  int chain[SIZE * SIZE][2];
  int reached[SIZE * SIZE][2];
  State color;
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      for (int ind = 0; ind < SIZE * SIZE; ++ind) {
        chain[ind][0] = -1;
        reached[ind][0] = -1;
        chain[ind][1] = -1;
        reached[ind][1] = -1;
      }
      if (board[i][j] == EMPTY && !managed[i][j]) {
        d_findReached(board, i, j, chain, reached);
        int chain_size = 0;
        while (chain[chain_size][0] >= 0 && chain[chain_size][0] < SIZE &&
               chain[chain_size][1] >= 0 && chain[chain_size][1] < SIZE) {
          managed[chain[chain_size][0]][chain[chain_size][1]] = true;
          ++chain_size;
        }
        color = board[reached[0][0]][reached[0][1]];
        int tmp = 0;
        while (reached[tmp][0] != -1) {
          if (board[reached[tmp][0]][reached[tmp][1]] != color) {
            break;
          }
          ++tmp;
        }
        if (reached[tmp][0] == -1) {
          if (color == BLACK) {
            results[0] += chain_size;
          } else {
            results[1] += chain_size;
          }
        }
      }
    }
  }

  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (board[i][j] == BLACK) {
        results[0] = results[0] + 1;
      } else if (board[i][j] == WHITE) {
        results[1] = results[1] + 1;
      }
    }
  }
}

__global__ void randomPlaysKernel(
    State *d_flattenedCubes,
    double *d_black_scores, // kazdy blok wpisze liczbe punktów czarnych po
                            // wszystkich symulacjach dla odpowiadającego diecka
    int *d_taken_black_stones, int *d_taken_white_stones,
    State *state_in_simulation) {
  __shared__ double
      localCounts[MAX_NUMBER_OF_THREADS]; // kazdy wątek zazanacza w swojej
                                          // komórce wynik symulacji dla
                                          // czarnych (1, 0.5 lub 0)
  int taken_stones[SIZE * SIZE][2];
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  curandState cs;
  curand_init(clock64() + tid, 0, 0, &cs);
  State board_for_random_play[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      board_for_random_play[i][j] =
          d_flattenedCubes[SIZE * SIZE * blockIdx.x + SIZE * i + j];
    }
  }
  State state = BLACK;
  if (*state_in_simulation == BLACK) {
    state = WHITE;
  }

  int num_of_tries = 0;
  int index = 0;
  int lost_black_stones = d_taken_black_stones[blockIdx.x];
  int lost_white_stones = d_taken_white_stones[blockIdx.x];
  int random_row, random_col;
  bool could_place_stone;

  while (index < NUM_OF_MOVEMENTS_IN_SIMULATION) {
    for (int i = 0; i < SIZE * SIZE; ++i) {
      taken_stones[i][0] = -1;
      taken_stones[i][1] = -1;
    }
    do {
      ++num_of_tries;
      random_row = curand(&cs) % SIZE;
      random_col = curand(&cs) % SIZE;

      could_place_stone = d_couldPlaceStone(board_for_random_play, random_row,
                                            random_col, state, taken_stones);
    } while (!could_place_stone && num_of_tries < 20);
    if (num_of_tries == 20) {
      break;
    }
    int num_of_taken_stones = 0;
    while (taken_stones[num_of_taken_stones][0] != -1) {
      ++num_of_taken_stones;
    }

    board_for_random_play[random_row][random_col] = state;
    if (state == BLACK) {
      lost_white_stones += num_of_taken_stones;
    } else {
      lost_black_stones += num_of_taken_stones;
    }
    if (state == BLACK) {
      state = WHITE;
    } else {
      state = BLACK;
    }
    int ind = 0;
    while (taken_stones[ind][0] != -1) {
      board_for_random_play[taken_stones[ind][0]][taken_stones[ind][1]] = EMPTY;
      ++ind;
    }
    num_of_tries = 0;
    ++index;
  }

  int results[2] = {0, 0};
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (board_for_random_play[i][j] == BLACK) {
        results[0] = results[0] + 1;
      } else if (board_for_random_play[i][j] == WHITE) {
        results[1] = results[1] + 1;
      }
    }
  }
  if (state == WHITE) {
    results[1] = results[1] + 1;
  }
  //__syncthreads();
  if ((results[0] + lost_white_stones) > (results[1] + lost_black_stones)) {
    localCounts[threadIdx.x] = 1.0;
  } else if ((results[0] + lost_white_stones) ==
             (results[1] + lost_black_stones)) {
    localCounts[threadIdx.x] = 0.5;
  }
  __syncthreads();
  if (threadIdx.x ==
      0) { // jeden główny wątek, który zapisze wuynik od całego bloku
    for (int i = 0; i < MAX_NUMBER_OF_THREADS; ++i) {
      d_black_scores[blockIdx.x] += localCounts[i];
    }
  }
}

void simulate(Node *n, State state) {
  cudaError_t cudaStatus;
  int totalSize = n->children.size() * SIZE * SIZE;
  State *h_flattenedCubes = new State[totalSize];
  flattenCube(n, h_flattenedCubes);
  State *d_flattenedCubes;
  cudaStatus = cudaMalloc(&d_flattenedCubes, totalSize * sizeof(State));
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMalloc (d_flattenedCubes) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }

  std::chrono::high_resolution_clock::time_point start, end;
  start = std::chrono::high_resolution_clock::now();
  cudaStatus = cudaMemcpy(d_flattenedCubes, h_flattenedCubes,
                          totalSize * sizeof(State), cudaMemcpyHostToDevice);
  end = std::chrono::high_resolution_clock::now();
  total_time_CPU_GPU_copying +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMemcpy (d_flattenedCubes) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  double *h_black_scores = new double[n->children.size()];
  double *d_black_scores; // kazde dziecko ma NUM_OF_THREADS
                          // symulacji, tu kazde zapisze liczba
                          // wygranych dla czarnych
  cudaStatus = cudaMalloc(&d_black_scores, n->children.size() * sizeof(double));
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMalloc (d_black_scores) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  cudaStatus =
      cudaMemset(d_black_scores, 0, n->children.size() * sizeof(double));
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMemset (d_black_scores) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  // state_in_simulation = state;
  State *d_state;
  cudaStatus = cudaMalloc(&d_state, sizeof(State));
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMalloc (d_state) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  start = std::chrono::high_resolution_clock::now();
  cudaStatus =
      cudaMemcpy(d_state, &state, sizeof(State), cudaMemcpyHostToDevice);
  end = std::chrono::high_resolution_clock::now();
  total_time_CPU_GPU_copying +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMemcpy (d_state) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  int *h_taken_white_stones = new int[n->children.size()];
  int *d_taken_white_stones;
  cudaStatus = cudaMalloc((void **)&d_taken_white_stones,
                          n->children.size() * sizeof(int));
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMalloc (d_taken_white_stones) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  int *h_taken_black_stones = new int[n->children.size()];
  int *d_taken_black_stones;
  cudaStatus = cudaMalloc((void **)&d_taken_black_stones,
                          n->children.size() * sizeof(int));
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMalloc (d_taken_black_stones) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  for (int i = 0; i < n->children.size(); ++i) {
    h_taken_black_stones[i] = n->children[i]->taken_black_stones;
    h_taken_white_stones[i] =
        n->children[i]->taken_white_stones; // juz zdobyte kamienie
  }
  start = std::chrono::high_resolution_clock::now();
  cudaStatus =
      cudaMemcpy(d_taken_white_stones, h_taken_white_stones,
                 n->children.size() * sizeof(int), cudaMemcpyHostToDevice);

  end = std::chrono::high_resolution_clock::now();
  total_time_CPU_GPU_copying +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMemcpy (d_taken_white_stones) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }
  start = std::chrono::high_resolution_clock::now();
  cudaStatus =
      cudaMemcpy(d_taken_black_stones, h_taken_black_stones,
                 n->children.size() * sizeof(int), cudaMemcpyHostToDevice);
  end = std::chrono::high_resolution_clock::now();
  total_time_CPU_GPU_copying +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  if (cudaStatus != cudaSuccess) {
    std::cout << "[ERROR] cudaMemcpy (d_taken_black_stones) failed: "
              << cudaGetErrorString(cudaStatus) << std::endl;
    exit(1);
  }

  randomPlaysKernel<<<n->children.size(), MAX_NUMBER_OF_THREADS>>>(
      d_flattenedCubes, d_black_scores, d_taken_black_stones,
      d_taken_white_stones, d_state);
  cudaDeviceSynchronize();
  start = std::chrono::high_resolution_clock::now();
  cudaMemcpy(h_black_scores, d_black_scores,
             n->children.size() * sizeof(double), cudaMemcpyDeviceToHost);
  end = std::chrono::high_resolution_clock::now();
  total_time_CPU_GPU_copying +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  for (int i = 0; i < n->children.size(); ++i) {
    n->number_of_simulations += MAX_NUMBER_OF_THREADS;
    n->children[i]->number_of_simulations += MAX_NUMBER_OF_THREADS;
    n->black_score = n->black_score + h_black_scores[i];
    n->children[i]->black_score =
        n->children[i]->black_score + h_black_scores[i];
  }
  delete[] h_flattenedCubes;
  delete[] h_black_scores;
  delete[] h_taken_black_stones;
  delete[] h_taken_white_stones;
  cudaFree(d_flattenedCubes);
  cudaFree(d_black_scores);
  cudaFree(d_taken_black_stones);
  cudaFree(d_taken_white_stones);
  cudaFree(d_state);
}

Node *findMaxUctChild(Node *parent, State state) {
  double maxUCT = -std::numeric_limits<double>::infinity();
  Node *maxUCTChild = parent->children[0];
  std::vector<Node *> topChildren = std::vector<Node *>();
  if (state == BLACK) {
    topChildren = parent->getMaxProfitForBlack();
  } else {
    topChildren = parent->getMaxProfitForWhite();
  }
  for (Node *child : topChildren) {
    double uctValue = calculateUct(child, state);
    if (uctValue > maxUCT) {
      maxUCT = uctValue;
      maxUCTChild = child;
    }
  }
  return maxUCTChild;
}

void backpropagate(Node *n) {
  Node *tmp = n;
  while (tmp->parent) {
    tmp->parent->black_score += tmp->black_score;
    tmp->parent->number_of_simulations += tmp->number_of_simulations;
    tmp = tmp->parent;
  }
}

Node *makeHumanMove(Node *parent, State state, int i, int j) {
  for (Node *child : parent->children) {
    if (child->board[i][j] == state) {
      return child;
    }
  }
  std::cout << "Helloo\n";
  return parent->children[0];
}

void showResults(Node *root_node, State actual_state) {
  std::cout << "Now we will see the results\n";
  std::cout << "Previous position for black:\n";
  printBoard(previousPositionForBlack);
  std::cout << "Previous position for white:\n";
  printBoard(previousPositionForWhite);
  if (actual_state == BLACK) {
    copyBoard(previousPositionForBlack, root_node->board);
  } else {
    copyBoard(previousPositionForWhite, root_node->board);
  }
  std::cout << "Main board\n";
  printBoard(root_node->board);
  auto main_results = computeTerritories(root_node->board);
  std::cout << "\nBlack territory: " << main_results.first << '\n';
  std::cout << "White territory: " << main_results.second << '\n';
  int lost_black_stones = total_taken_black;
  int lost_white_stones = total_taken_white;
  std::cout << "Lost black stones: " << total_taken_black << '\n';
  std::cout << "Lost white stones: " << total_taken_white << '\n';
  std::cout << "-------------------------------------\n";
  if ((main_results.first + lost_white_stones) >
      (main_results.second + lost_black_stones)) {
    std::cout << "BLACK won\n";
  } else if ((main_results.first + lost_white_stones) <
             (main_results.second + lost_black_stones)) {
    std::cout << "WHITE won\n";
  } else {
    std::cout << "DRAW\n";
  }
  std::cout << "-------------------------------------\n";
}

void play(Node *root_node, State actual_state, bool isHumanVsComp,
          State humanState) {
  Node *actual_node;
  State whoose_move = actual_state;
  int max_depth_ind = 0;
  int mov_ind = 0;
  int row_by_user, col_by_user;
  std::cout << "\nStart board: \n";
  printBoard(root_node->board);
  std::chrono::high_resolution_clock::time_point start, end;
  while (mov_ind < MOVEMENTS) {
    max_depth_ind = 0;
    while (max_depth_ind < MAX_DEPTH) {
      actual_node = root_node;
      int local_depth = 0;
      while (actual_node->children.size() != 0) {
        if (local_depth % 2 == 0) {
          whoose_move = actual_state;
        } else {
          whoose_move = changeState(actual_state);
        }
        start = std::chrono::high_resolution_clock::now();
        actual_node = findMaxUctChild(actual_node, whoose_move); // select
        end = std::chrono::high_resolution_clock::now();
        ++selection_moves;
        total_time_selection +=
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        ++local_depth;
      }
      if (local_depth % 2 == 0) {
        whoose_move = actual_state;
      } else {
        whoose_move = changeState(actual_state);
      }
      start = std::chrono::high_resolution_clock::now();
      expand(actual_node, whoose_move);
      end = std::chrono::high_resolution_clock::now();
      ++expansion_moves;
      total_time_expansion +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      start = std::chrono::high_resolution_clock::now();
      simulate(actual_node, whoose_move);
      end = std::chrono::high_resolution_clock::now();
      ++simulation_moves;
      total_time_simulation +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      start = std::chrono::high_resolution_clock::now();
      backpropagate(actual_node);
      end = std::chrono::high_resolution_clock::now();
      ++backpropagation_moves;
      total_time_backpropagation +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      ++max_depth_ind;
    }

    if (isHumanVsComp && actual_state == humanState) {
      std::cout << "Your move:\n";
      std::cin >> row_by_user >> col_by_user;
      if (row_by_user == -1)
        return;
      root_node =
          makeHumanMove(root_node, actual_state, row_by_user, col_by_user);
    } else {
      root_node = findMaxUctChild(root_node,
                                  actual_state);
    }

    std::cout << "\nNr: " << mov_ind << '\n';
    printBoard(root_node->board);
    total_taken_black = root_node->taken_black_stones;
    total_taken_white = root_node->taken_white_stones;
    std::cout << "Lost black stones: " << root_node->taken_black_stones << '\n';
    std::cout << "Lost white stones: " << root_node->taken_white_stones << '\n';
    if (actual_state == BLACK) { // przekazujemy prawo ruchu innemy graczowi
      actual_state = WHITE;
      copyBoard(root_node->board, previousPositionForWhite);
    } else {
      actual_state = BLACK;
      copyBoard(root_node->board, previousPositionForBlack);
    }
    ++mov_ind;
  }
}

void emptyBoard(State actual_board[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      actual_board[i][j] = EMPTY;
    }
  }
}

void preProcessing(Node *root_node, State &actual_state,
                   State actual_board[SIZE][SIZE], bool &is_black,
                   bool &isHumanVsComp, State &humanState, bool is_load_board,
                   std::string &filename) {
  std::srand(std::time(0));
  emptyBoard(actual_board);
  createNeighbours();

  int num_of_o = 0;
  int num_of_x = 0;
  if (is_load_board) {
    std::ifstream file(filename);
    if (file.is_open()) {
      std::string line;
      int i = 0;
      while (getline(file, line)) {
        int j = 0;
        for (char c : line) {
          if (c == 'o' || c == 'O') {
            actual_board[i][j] = WHITE;
            ++num_of_o;
          } else if (c == 'x' || c == 'X') {
            actual_board[i][j] = BLACK;
            ++num_of_x;
          }
          ++j;
        }
        ++i;
      }

      file.close();
    } else {
      std::cerr << "Unable to open file: " << filename << std::endl;
      exit(1);
    }
  }
  if (num_of_o != num_of_x) {
    actual_state = WHITE;
  }

  copyBoard(actual_board, root_node->board);
  std::cout << "Input board:\n";
  printBoard(root_node->board);
  if (actual_state == BLACK) {
    copyBoard(actual_board, previousPositionForBlack);
    emptyBoard(previousPositionForWhite);
  } else {
    copyBoard(actual_board, previousPositionForWhite);
    emptyBoard(previousPositionForBlack);
  }

  int tmp;
  std::cout
      << "Select mode:\n 1 - copmuter vs computer\n 2 - human vs computer\n";
  std::cin >> tmp;
  if (tmp == 2) {
    isHumanVsComp = true;
  }
  if (isHumanVsComp) {
    std::cout << "Select color:\n 1 - black\n 2 - white\n";
    std::cin >> tmp;
    if (tmp == 2) {
      humanState = WHITE;
    }
  }
}

void showTime() {
  double average_time_selection =
      total_time_selection.count() / (double)selection_moves;
  double average_time_expansion =
      total_time_expansion.count() / (double)expansion_moves;
  double average_time_simulation =
      total_time_simulation.count() / (double)simulation_moves;
  double average_time_backpropagation =
      total_time_backpropagation.count() / (double)backpropagation_moves;
  std::cout << "Selection average time: " << average_time_selection << " \n";
  std::cout << "Expansion average time: " << average_time_expansion << " \n";
  std::cout << "Simulation average time: " << average_time_simulation << " \n";
  std::cout << "Backpropagation average time: " << average_time_backpropagation
            << " \n";
  std::cout << "CPU-GPU data copying: " << total_time_CPU_GPU_copying.count()
            << '\n';
}

int main(int argc, char **argv) {
  State actual_state = BLACK;
  bool is_black = true;
  State actual_board[SIZE][SIZE];
  Node *root_node = new Node;
  bool isHumanVsComp = false;
  bool is_load_board = false;
  std::string filename = "";
  if (argc >= 2) {
    is_load_board = true;
    filename = argv[1];
  }
  State humanState = BLACK;
  preProcessing(root_node, actual_state, actual_board, is_black, isHumanVsComp,
                humanState, is_load_board, filename);
  play(root_node, actual_state, isHumanVsComp, humanState);
  showResults(root_node, actual_state);
  showTime();
  delete root_node;
  return 0;
}
