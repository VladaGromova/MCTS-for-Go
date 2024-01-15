#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <math.h>
#include <string>
#include <utility>
#include <vector>

#define C sqrt(2)
#define SIZE 9
#define NUM_OF_MOVEMENTS_IN_SIMULATION 10
#define MAX_DEPTH 15 // tyle razy wykonamy te 3 etapy
#define MOVEMENTS 60

enum State { EMPTY, BLACK, WHITE };

std::vector<std::pair<int, int>> NEIGHBOURS[SIZE][SIZE];
State previousPositionForBlack[SIZE][SIZE];
State previousPositionForWhite[SIZE][SIZE];

typedef struct Node {
  State board[SIZE][SIZE];
  std::vector<Node *> children = std::vector<Node *>();
  Node *parent = NULL;
  int taken_black_stones = 0;
  int taken_white_stones = 0;
  unsigned int number_of_simulations = 0;
  double black_score = 0.0;
  double uct;
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
    return std::numeric_limits<double>::
        infinity(); // Return infinity if the child has not been explored yet.
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
      if (prev[i][j] != actual[i][j])
        return false;
    }
  }
  return true;
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
      } else if (board[fn.first][fn.second] != color &&
                 std::find(reached.begin(), reached.end(), fn) ==
                     reached.end()) {
        reached.push_back(fn);
      }
    }
  }
  return std::make_pair(chain, reached);
}

std::pair<bool, std::vector<std::pair<int, int>>>
couldPlaceStone(State board[SIZE][SIZE], int row, int col, State state) {
  std::vector<std::pair<int, int>> taken_stones =
      std::vector<std::pair<int, int>>();
  if (board[row][col] != EMPTY)
    return std::make_pair(false, taken_stones);
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
          // board[potential_captured[ind].first][potential_captured[ind].second]
          // = EMPTY; // NEW
          taken_stones.push_back(std::make_pair(
              potential_captured[ind].first, potential_captured[ind].second));
        }
        /* taken_stones.push_back(std::make_pair(
             stones_reached_by_mine[i].first,
           stones_reached_by_mine[i].second));*/
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

void createChildren(Node *n, int i, int j, State state) {
  // Node childNode;
  auto taken_stones = couldPlaceStone(n->board, i, j, state).second;
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

State randomPlays(Node *n, State state) {
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
  int lost_white_stones = n->taken_white_stones;
  int lost_black_stones = n->taken_black_stones;
  while (index < NUM_OF_MOVEMENTS_IN_SIMULATION) {
    ++index;
    do {
      ++num_of_tries;
      generateRandomCell(random_row, random_col);
      result =
          couldPlaceStone(board_for_random_play, random_row, random_col, state);
    } while (!result.first && num_of_tries < 20);
    if (num_of_tries == 20)
      break;
    std::vector<std::pair<int, int>> taken_stones = result.second;
    board_for_random_play[random_row][random_col] = state;
    if (state == BLACK) {
      lost_white_stones += taken_stones.size();
      copyBoard(board_for_random_play, previousPositionForWhite);
    } else {
      lost_black_stones += taken_stones.size();
      copyBoard(board_for_random_play, previousPositionForBlack);
    }
    // change state from black to white and vice versa
    if (state == BLACK) {
      state = WHITE;
    } else {
      state = BLACK;
    }
    for (int i = 0; i < taken_stones.size(); ++i) { // taking oppoents' stones
      board_for_random_play[taken_stones[i].first][taken_stones[i].second] =
          EMPTY;
    }
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
            state); // i j - tu ustawimy kamien i takie dziecko dodamy do n
      }
    }
  }
}

void simulate(Node *n, State state) {
  State winner;
  for (int i = 0; i < n->children.size(); ++i) {
    if (state == BLACK) {
      copyBoard(n->board, previousPositionForBlack);
    } else {
      copyBoard(n->board, previousPositionForWhite);
    }
    winner = randomPlays(n->children[i], state);
    n->number_of_simulations = n->number_of_simulations + 1;
    n->children[i]->number_of_simulations =
        n->children[i]->number_of_simulations + 1;
    if (winner == BLACK) {
      n->black_score = n->black_score + 1.0;
      n->children[i]->black_score = n->children[i]->black_score + 1.0;
    } else if (winner == EMPTY) {
      n->black_score = n->black_score + 0.5;
      n->children[i]->black_score = n->children[i]->black_score + 0.5;
    }
  }
}



Node *findMaxUctChild(Node *parent, State state) {
  double maxUCT = -std::numeric_limits<double>::infinity();
  Node *maxUCTChild = nullptr;
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

void testing() {
  createNeighbours();
  State testing_board[SIZE][SIZE] = {{BLACK, WHITE, WHITE, EMPTY, EMPTY},
                                     {EMPTY, BLACK, BLACK, WHITE, EMPTY},
                                     {EMPTY, EMPTY, EMPTY, EMPTY, EMPTY},
                                     {EMPTY, EMPTY, EMPTY, EMPTY, EMPTY},
                                     {BLACK, EMPTY, EMPTY, EMPTY, WHITE}};
  auto res = couldPlaceStone(testing_board, 0, 3, BLACK);
  int i = 0;
}

Node *makeHumanMove(Node *parent, State state, int i, int j) {
  for (Node *child : parent->children) {
    if (child->board[i][j] == state) {
      return child;
    }
  }
}

int main(int argc, char **argv) {
  // testing();
  int tmp;
  State actual_state = BLACK;
  std::srand(std::time(0));
  createNeighbours();
  bool is_black = true;
  State actual_board[SIZE][SIZE];
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      actual_board[i][j] = EMPTY;
    }
  }
  std::cout << "Do you want to load board? 1 - yes, 2 - no\n";
  std::cin >> tmp;
  
  
  if (tmp == 1) {
    std::cin.ignore();
    std::string input;
    int num_of_x = 0;
    int num_of_o = 0;
    for (int i = 0; i < SIZE; ++i) {
      std::getline(std::cin, input);
      int char_nr = 0;
      for (char c : input) {
        if (c == 'X' || c == 'x') {
          ++num_of_x;
          actual_board[i][char_nr] = BLACK;
        } else if (c == 'O' || c == 'o') {
          ++num_of_o;
          actual_board[i][char_nr] = WHITE;
        }
        ++char_nr;
      }
    }
    if(num_of_x != num_of_o){
      actual_state = WHITE;
    }
  }

  Node MCTS_head;
  copyBoard(actual_board, MCTS_head.board);
  copyBoard(actual_board, previousPositionForBlack);
  copyBoard(actual_board, previousPositionForWhite);
  Node *actual_node;
  Node *root_node;
  root_node = &MCTS_head;   // aktualny stan pozycji
  actual_node = &MCTS_head; // will be changed later (selection)

  int num_of_test_movements = 0;
  int num_of_tree_searches = 0;
  int lost_white_stones = 0;
  int lost_black_stones = 0;
  State whoose_move = actual_state;
  int max_depth_ind = 0;
  int mov_ind = 0;

  bool isHumanVsComp = false;
  State humanState = BLACK;
  int row_by_user, col_by_user;
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

  
    std::cout << "\nStart board: \n";
  printBoard(root_node);
  while (mov_ind < MOVEMENTS) {
    max_depth_ind = 0;
    while (max_depth_ind < MAX_DEPTH) { 
      actual_node = root_node;
      int local_depth = 0;
      while (actual_node->children.size() != 0) {
        if (local_depth % 2 == 0) {
          whoose_move = actual_state;
        } else {
          if (actual_state ==
              BLACK) {  // przekazujemy prawo ruchu innemy graczowi
            whoose_move = WHITE;
          } else {
            whoose_move = BLACK;
          }
        }
        actual_node = findMaxUctChild(actual_node, whoose_move); // select
        ++local_depth;
      }
      if (local_depth % 2 == 0) {
        whoose_move = actual_state;
      } else {
        if (actual_state == BLACK) { // przekazujemy prawo ruchu innemy graczowi
          whoose_move = WHITE;
        } else {
          whoose_move = BLACK;
        }
      }
      expand(actual_node, whoose_move);
      simulate(actual_node, whoose_move);
      backpropagate(actual_node);
      ++max_depth_ind;
    }

    if (isHumanVsComp && actual_state == humanState) {
      std::cout << "Your move:\n";
      std::cin >> row_by_user >> col_by_user;
      root_node =
          makeHumanMove(root_node, actual_state, row_by_user, col_by_user);
    } else {
      root_node = findMaxUctChild(root_node, actual_state); // na pewno bo tu robimy ruch
    }

    std::cout << "\nNr: " << mov_ind << '\n';
    printBoard(root_node);
    std::cout << "Lost black stones: " << root_node->taken_black_stones << '\n';
    std::cout << "Lost white stones: " << root_node->taken_white_stones << '\n';
    if (actual_state == BLACK) { // przekazujemy prawo ruchu innemy graczowi
      actual_state = WHITE;
    } else {
      actual_state = BLACK;
    }
    ++mov_ind;
  }

  auto main_results = computeTerritories(root_node->board);
  std::cout << "\nBlack territory: " << main_results.first << '\n';
  std::cout << "White territory: " << main_results.second << '\n';
  lost_black_stones = root_node->taken_black_stones;
  lost_white_stones = root_node->taken_white_stones;
  std::cout << "Lost black stones: " << lost_black_stones << '\n';
  std::cout << "Lost white stones: " << lost_white_stones << '\n';
  if ((main_results.first + lost_white_stones) >
      (main_results.second + lost_black_stones)) {
    std::cout << "BLACK won\n";
  } else if ((main_results.first + lost_white_stones) <
             (main_results.second + lost_black_stones)) {
    std::cout << "WHITE won\n";
  } else {
    std::cout << "DRAW\n";
  }
  return 0;
}
