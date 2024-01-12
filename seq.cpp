#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <utility>
#include <vector>

#define SIZE 9
#define C sqrt(2)

enum State { EMPTY, BLACK, WHITE };

std::vector<std::pair<int, int>> NEIGHBOURS[SIZE][SIZE];
State previousPositionForBlack[SIZE][SIZE];
State previousPositionForWhite[SIZE][SIZE];

typedef struct Node {
  State board[SIZE][SIZE];
  std::vector<Node> children = std::vector<Node>();
  int lost_white_stones = 0;
  int lost_black_stones = 0;
  Node *parent = NULL;
  unsigned int number_of_simulations = 0;
  double black_score = 0.0;
} Node;

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
      std::cout << n->board[i][j] << '\t';
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

void copyBoard(Node *parent, Node *child) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      child->board[i][j] = parent->board[i][j];
    }
  }
}

void createChildren(Node *n, int i, int j, State state) {
  Node childNode;
  copyBoard(n, &childNode);
  childNode.board[i][j] = state;
  n->children.push_back(childNode);
  childNode.parent = n;
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
findReached(Node *node, int i, int j) {
  State color = node->board[i][j];
  std::vector<std::pair<int, int>> chain, reached;
  std::vector<std::pair<int, int>> frontier = {std::make_pair(i, j)};
  while (!frontier.empty()) {
    auto current_fc = frontier.back();
    frontier.pop_back();
    chain.push_back(current_fc);

    for (auto fn : NEIGHBOURS[i][j]) {
      if (node->board[i][j] == color &&
          std::find(chain.begin(), chain.end(), fn) == chain.end()) {
        frontier.push_back(fn);
      } else if (node->board[i][j] != color) {
        reached.push_back(fn);
      }
    }
  }
  return std::make_pair(chain, reached);
}

std::pair<bool, std::vector<std::pair<int, int>>>
couldPlaceStone(Node *n, int row, int col, State state) {
  std::vector<std::pair<int, int>> taken_stones =
      std::vector<std::pair<int, int>>();
  if (n->board[row][col] != EMPTY)
    return std::make_pair(false, taken_stones);
  n->board[row][col] = state;
  State alternativeState;
  if (state == BLACK) {
    alternativeState = WHITE;
  } else {
    alternativeState = BLACK;
  }
  std::vector<std::pair<int, int>> stones_reached_by_mine =
      findReached(n, row, col).second;
  // tu sprawdzamy czy cos zabieramy u przeciwnika
  for (int i = 0; i < stones_reached_by_mine.size(); ++i) {
    if (n->board[stones_reached_by_mine[i].first]
                [stones_reached_by_mine[i].second] !=
        EMPTY) { // dla kazego kamienia z sąsiadujących
      auto reached_by_opponent = findReached(n, stones_reached_by_mine[i].first,
                                             stones_reached_by_mine[i].second)
                                     .second;
      int j = 0;
      for (j = 0; j < reached_by_opponent.size();
           ++j) { // dla kazdego wezla co jest obok kamienia przeciwnika
        if (n->board[reached_by_opponent[j].first]
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
      if (n->board[stones_reached_by_mine[j].first][stones_reached_by_mine[j].second] != alternativeState) {
        break;
      }
    }
    if(j== stones_reached_by_mine.size()){ // only opponent stones surround me, suicide
      return std::make_pair(false, taken_stones);
    }
  }
  // end
  // sprawdzamy czy nie ko
  // TO DO
  //
  n->board[row][col] = EMPTY;
  return std::make_pair(true, taken_stones);
}

void simulate(Node *n, State state) {
  bool in_process = true;
  int random_row, random_col;
  std::pair<bool, std::vector<std::pair<int, int>>> result;
  while (in_process) {
    do {
      generateRandomCell(random_row, random_col);
      result = couldPlaceStone(n, random_row, random_col, state);
    } while (!result.first);
    std::vector<std::pair<int, int>> taken_stones = result.second;
    n->board[random_row][random_col] = state;
    if(state == BLACK){
      n->lost_white_stones += taken_stones.size();
    } else {
      n->lost_black_stones += taken_stones.size();
    }
    for(int i=0; i<taken_stones.size(); ++i){
      n->board[taken_stones[i].first][taken_stones[i].second] = EMPTY;
    }
    // TO DO: save previous position somewhere

    // change state from black to white and vice versa
    //
  }
}

void expand(Node *n, State state) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      if (couldPlaceStone(n, i, j, state).first) {
        createChildren(
            n, i, j,
            state); // i j - tu ustawimy kamien i takie dziecko dodamy do n
      }
    }
  }
  for (int i = 0; i < n->children.size(); ++i) {
    simulate(&n->children[i], state);
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
  createRootNode(&MCTS_head, actual_board);
  printBoard(&MCTS_head);
  // while....
  Node actual_node;
  actual_node = MCTS_head; // will be changed later (selection)
  expand(&actual_node, State::BLACK);

  std::cout << "test\n";
  return 0;
}
