#include <math.h>

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <utility>
#include <vector>

#define SIZE 5
#define C sqrt(2)

enum State { EMPTY, BLACK, WHITE };

std::vector<std::pair<int, int>> NEIGHBOURS[SIZE][SIZE];
State previousPositionForBlack[SIZE][SIZE];
State previousPositionForWhite[SIZE][SIZE];

typedef struct Node {
  State board[SIZE][SIZE];
  std::vector<Node *> children = std::vector<Node *>();
  Node *parent = NULL;
  unsigned int number_of_simulations = 0;
  double black_score = 0.0;
  double uct;
  ~Node() {
    for (Node *child : children) {
      delete child;
    }
  }
} Node;

double calculateUct(Node *n, State state) {
  if (n->number_of_simulations == 0) {
    return std::numeric_limits<double>::
        infinity();  // Return infinity if the child has not been explored yet.
  }

  // UCT formula
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

void createRootNode(Node *n, State b[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      n->board[i][j] = b[i][j];
    }
  }
}

void printBoard(Node *n) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (n->board[i][j] == EMPTY) {
        std::cout << ".\t";
      } else if (n->board[i][j] == BLACK) {
        std::cout << "X\t";
      } else {
        std::cout << "O\t";
      }
    }
    std::cout << '\n';
  }
  std::cout << '\n';
}

bool isKo(State prev[SIZE][SIZE], State actual[SIZE][SIZE]) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (prev[i][j] != actual[i][j]) return false;
    }
  }
  return true;
}

void createChildren(Node *n, int i, int j, State state) {
  //Node childNode;
  Node *childNode = new Node(); 
  copyBoard(n->board, childNode->board);
  childNode->board[i][j] = state;
  n->children.push_back(childNode);
  childNode->parent = n;
}

void generateRandomCell(int &randomRow, int &randomCol) {
  randomRow = std::rand() % SIZE;
  randomCol = std::rand() % SIZE;
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
      } else if (board[fn.first][fn.second] != color && std::find(reached.begin(), reached.end(), fn) == reached.end()) {
        reached.push_back(fn);
      }
    }
  }
  return std::make_pair(chain, reached);
}

std::pair<bool, std::vector<std::pair<int, int>>> couldPlaceStone(
    State board[SIZE][SIZE],
                                                                  int row,
                                                                  int col,
                                                                  State state) {
  std::vector<std::pair<int, int>> taken_stones =
      std::vector<std::pair<int, int>>();
  if (board[row][col] != EMPTY) return std::make_pair(false, taken_stones);
  board[row][col] = state;
  State alternativeState;
  if (state == BLACK) {
    alternativeState = WHITE;
  } else {
    alternativeState = BLACK;
  }
  std::vector<std::pair<int, int>> stones_reached_by_mine =
      findReached(board, row, col).second;
  // tu sprawdzamy czy cos zabieramy u przeciwnika
  for (int i = 0; i < stones_reached_by_mine.size(); ++i) {
    if (board[stones_reached_by_mine[i].first]
                [stones_reached_by_mine[i].second] !=
        EMPTY) {  // dla kazego kamienia z sąsiadujących
      auto reached_by_opponent = findReached(board, stones_reached_by_mine[i].first,
                                             stones_reached_by_mine[i].second)
                                     .second;
      int j = 0;
      for (j = 0; j < reached_by_opponent.size();
           ++j) {  // dla kazdego wezla co jest obok kamienia przeciwnika
        if (board[reached_by_opponent[j].first]
                    [reached_by_opponent[j].second] != state) {
          break;
        }
      }
      if (j == reached_by_opponent.size()) {
        taken_stones.push_back(std::make_pair(
            stones_reached_by_mine[i].first, stones_reached_by_mine[i].second));
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
                 .size()) {  // only opponent stones surround me, suicide
      return std::make_pair(false, taken_stones);
    }
  }
  // end
  // sprawdzamy czy nie ko
  if (state == BLACK) {
    isKo(previousPositionForBlack, board);
  } else {
    isKo(previousPositionForWhite, board);
  }
  
  board[row][col] = EMPTY;
  return std::make_pair(true, taken_stones);
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
        // jesli w reached wszystkie kamienie sa jednego koloru to to jest terytorium tego koloru
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

State simulate(Node *n, State state) {
  State board_for_random_play[SIZE][SIZE];
  copyBoard(n->board, board_for_random_play);
  int random_row, random_col;
  std::pair<bool, std::vector<std::pair<int, int>>> result;
  if (state == BLACK) {
    state = WHITE;
    copyBoard(n->board, previousPositionForWhite);
  } else {
    state = BLACK;
    copyBoard(n->board, previousPositionForBlack);
  }
  int num_of_tries = 0;
  int index = 0;
  int lost_white_stones = 0;
  int lost_black_stones = 0;
  while (index<20) {
    ++index;
    do {
      ++num_of_tries;
      generateRandomCell(random_row, random_col);
      result = couldPlaceStone(board_for_random_play, random_row, random_col, state);
    } while (!result.first && num_of_tries < 20);
    if (num_of_tries == 20) break;
    std::vector<std::pair<int, int>> taken_stones = result.second;
    board_for_random_play[random_row][random_col] = state;
    if (state == BLACK) {
      lost_white_stones += taken_stones.size();
      copyBoard(board_for_random_play, previousPositionForWhite);
    } else {
      lost_black_stones += taken_stones.size();
      copyBoard(board_for_random_play, previousPositionForBlack);
    }
    for (int i = 0; i < taken_stones.size(); ++i) { // taking oppoents' stones
      board_for_random_play[taken_stones[i].first][taken_stones[i].second] =
          EMPTY;
    }
    // change state from black to white and vice versa
    if (state == BLACK) {
      state = WHITE;
    } else {
      state = BLACK;
    }
    //
    // !!! create condition for in_process (do poki są EMPTY ??)
    num_of_tries = 0;
  }
  // count scores and territories
  auto results = computeTerritories(board_for_random_play);
  if ((results.first + lost_white_stones) >
      (results.second + lost_black_stones)) {
    return BLACK;
  } else if ((results.first + lost_white_stones) <
             (results.second + lost_black_stones)) {
    return WHITE;
  } else {
    return EMPTY;
  }
}

void expand(Node *n, State state) {
  if (state == BLACK) {
    copyBoard(n->board, previousPositionForBlack);
  } else {
    copyBoard(n->board, previousPositionForWhite);
  }
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (couldPlaceStone(n->board, i, j, state).first) {
        createChildren(
            n, i, j,
            state);  // i j - tu ustawimy kamien i takie dziecko dodamy do n
      }
    }
  }
  State winner;
  for (int i = 0; i < n->children.size(); ++i) {
    if (state == BLACK) {
      copyBoard(n->board, previousPositionForBlack);
    } else {
      copyBoard(n->board, previousPositionForWhite);
    }
    winner = simulate(n->children[i], state);
    n->number_of_simulations = n->number_of_simulations + 1;
    if (winner == BLACK) {
      n->black_score = n->black_score + 1.0;
    } else if (winner == EMPTY) {
      n->black_score = n->black_score + 0.5;
    }
  }
}

Node *findMaxUctChild(Node *parent, State state) {
  double maxUCT = -std::numeric_limits<double>::infinity();
  Node *maxUCTChild = nullptr;

  for (Node *child : parent->children) {
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
  while (tmp->parent){
    tmp->parent->black_score += tmp->black_score;
    tmp->parent->number_of_simulations += tmp->number_of_simulations;
    tmp = tmp->parent;
  }
}

int main(int argc, char **argv) {
  std::srand(std::time(0));
  createNeighbours();
  bool is_black = true;
  State actual_board[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      actual_board[i][j] = EMPTY;
    }
  }
  
  Node MCTS_head;
  copyBoard(actual_board, MCTS_head.board);
  copyBoard(actual_board, previousPositionForBlack);
  copyBoard(actual_board, previousPositionForWhite);
  Node *actual_node;
  Node *root_node;
  root_node = &MCTS_head; // aktualny stan pozycji
  actual_node = &MCTS_head;  // will be changed later (selection)
  State actual_state = BLACK;

  int num_of_test_movements = 0;
  int num_of_tree_searches = 0;
  while (num_of_test_movements < 13) {
    ++num_of_test_movements;
    num_of_tree_searches = 0;
    while (num_of_tree_searches < 5) {
      ++num_of_tree_searches;
      actual_node = root_node;
      while (!actual_node->children.empty()) {  // select
        actual_node = findMaxUctChild(actual_node, actual_state);
      }
      expand(actual_node, actual_state);
      backpropagate(actual_node);
    }
    root_node = findMaxUctChild(root_node, actual_state); // making move
    printBoard(root_node);
    // !!!!!!!! change previous position somewhere
    if (actual_state == BLACK) {
      actual_state = WHITE;
    } else {
      actual_state = BLACK;
    }
  }

  std::cout << "test\n";
  return 0;
}

