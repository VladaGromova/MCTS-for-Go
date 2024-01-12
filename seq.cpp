#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#define SIZE 9
#define C sqrt(2)

enum State { EMPTY, BLACK, WHITE };

typedef struct Node {
  State board[SIZE][SIZE];
  std::vector<Node> children = std::vector<Node>();
  Node* parent = NULL;
  unsigned int number_of_simulations = 0;
  double black_score = 0.0;
} Node;

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

void copyBoard(Node *parent, Node *child) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      child->board[i][j] = parent->board[i][j];
    }
  }
}

void createChildren(Node *n, int i, int j, bool is_black) {
  Node childNode;
  copyBoard(n, &childNode);
  if (is_black) {
    childNode.board[i][j] = BLACK;
  } else {
    childNode.board[i][j] = WHITE;
  }
  n->children.push_back(childNode);
  childNode.parent = n;
}

void generateRandomCell(int& randomRow, int& randomCol) {
    randomRow = std::rand() % SIZE;
    randomCol = std::rand() % SIZE;
}

void simulate(Node *n, bool is_black){
  // petla do wygranej - TO DO 
  int random_row, random_col;
  do {
    generateRandomCell(random_row, random_col);
  }while (n->board[random_row][random_col] != EMPTY);
  if(is_black){ 
    n->board[random_row][random_col] = BLACK;
  } else {
    n->board[random_row][random_col] = WHITE;
  }
  is_black = !is_black;
  //
}

void expand(Node *n, bool is_black) {
  for (int i = 0; i < SIZE; ++i) {
    for (int j = 0; j < SIZE; ++j) {
      // if no suicide and no ko-rule
      // but now i'm just checking if there is no another stone
      if (n->board[i][j] == EMPTY) {
        createChildren(n, i, j, is_black); // i j - tu ustawimy kamien i takie dziecko dodamy do n
      }
    }
  }
  for(int i=0; i<n->children.size(); ++i){
    simulate(&n->children[i], is_black);
  }
}

int main(int argc, char **argv) {
  std::srand(std::time(0));
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
  expand(&actual_node, is_black);

  std::cout << "test\n";
  return 0;
}
