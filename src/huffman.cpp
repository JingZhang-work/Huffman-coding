#include "bits.h"
#include "treenode.h"
#include "huffman.h"
#include "map.h"
#include "vector.h"
#include "priorityqueue.h"
#include "strlib.h"
#include "testing/SimpleTest.h"
using namespace std;

bool isConnector(EncodingTreeNode *connector) {
    return connector->zero != nullptr && connector->one != nullptr;
}

bool isLeaf(EncodingTreeNode* node) {
    return node->zero == nullptr && node->one == nullptr;
}

/**
 * Given a Queue<Bit> containing the compressed message bits and the encoding tree
 * used to encode those bits, decode the bits back to the original message text.
 *
 * You can assume that tree is a well-formed non-empty encoding tree and
 * bits queue contains a valid sequence of encoded bits.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 */
string decodeText(EncodingTreeNode* tree, Queue<Bit>& messageBits) {
    // there is nothing to decode.
    if (messageBits.isEmpty()) {
        return "";
    }

    string decodedMessage = "";

    EncodingTreeNode *cur = tree;

    // there is bit which needs to be decode.
    while (!messageBits.isEmpty()) {
        Bit direction = messageBits.dequeue();

        // follow down the tree
        if (direction == 0) {
            cur = cur->zero;
        } else {
            cur = cur->one;
        }

        // get the char on the leaf node and add it to the message
        if (isLeaf(cur)) {
            decodedMessage += cur->ch;
            cur = tree;
        }
    }

    return decodedMessage;
}

/**
 * Reconstruct encoding tree from flattened form Queue<Bit> and Queue<char>.
 *
 * You can assume that the input Queues are well-formed and represent
 * a valid encoding tree.
 */
EncodingTreeNode* unflattenTree(Queue<Bit>& treeBits, Queue<char>& treeLeaves) {
    // there is nothing need to be unflatten
    if (treeBits.isEmpty() || treeLeaves.isEmpty()) {
        return nullptr;
    }

    // the node is a leaf, put the char
    if (treeBits.dequeue() == 0) {
        EncodingTreeNode* leafNode = new EncodingTreeNode(treeLeaves.dequeue());
        return leafNode;
    } else {// the node is a connector, make it and begin recursion
        EncodingTreeNode* connector = new EncodingTreeNode(nullptr, nullptr);
        connector->zero = unflattenTree(treeBits, treeLeaves);
        connector->one = unflattenTree(treeBits,treeLeaves);
        return connector;
    }
}

/**
 * Decompress the given EncodedData and return the original text.
 *
 * You can assume the input data is well-formed and was created by a correct
 * implementation of compress.
 *
 * Your implementation may change the data parameter however you like. There
 * are no requirements about what it should look like after this function
 * returns.
 */
string decompress(EncodedData& data) {
    /* TODO: Fill in the rest of this function. */
    EncodingTreeNode *root = unflattenTree(data.treeBits, data.treeLeaves);
    string decodedmessage = decodeText(root, data.messageBits);
    deallocateTree(root);
    return decodedmessage;
}

Map<char, int> getCharsFrequency(string text) {
    Map<char, int> charsFrequency;
    for (char character : text) {
        charsFrequency[character] += 1;
    }
    return charsFrequency;
}

/**
 * Constructs an optimal Huffman coding tree for the given text, using
 * the algorithm described in lecture.
 *
 * Reports an error if the input text does not contain at least
 * two distinct characters.
 *
 * When assembling larger trees out of smaller ones, make sure to set the first
 * tree dequeued from the queue to be the zero subtree of the new tree and the
 * second tree as the one subtree.
 */
EncodingTreeNode* buildHuffmanTree(string text) {

    // build frequency map
    Map<char, int> charsFrequency = getCharsFrequency(text);
    if (charsFrequency.size() < 2) {
        error("At least two distinct chars.");
    }

    // build priority queue
    PriorityQueue<EncodingTreeNode*> charsTree;
    for (char character : charsFrequency.keys()) {
        EncodingTreeNode *node = new EncodingTreeNode(character); // initialize a node
        charsTree.enqueue(node, charsFrequency[character]); // enqueue the node according to the frequency on the map
    }

    // combine two nodes
    while (charsTree.size() > 1) {
        int zeroFre = charsTree.peekPriority();
        EncodingTreeNode *zero = charsTree.dequeue(); // get the first one

        int oneFre = charsTree.peekPriority();
        EncodingTreeNode *one = charsTree.dequeue(); // get the second one

        EncodingTreeNode *root = new EncodingTreeNode(zero, one); // create a new one by combining them
        int sumFre = oneFre + zeroFre;
        charsTree.enqueue(root, sumFre); //enqueue to the tree
    }

    return charsTree.dequeue();
}

void getCharCode(Map<char, Vector<Bit>>& charAndCode, EncodingTreeNode *tree, Vector<Bit>& code) {
    // if there is nothing in the tree, just return it
    if (tree == nullptr) {
        return;
    }

    // if it is a leaf, get the path and add to the map
    if (isLeaf(tree)) {
        charAndCode[tree->ch] = code;
        return;

    } else { // if it is a connector
        // go to the left
        code.add(0);
        // tree = tree->zero;
        getCharCode(charAndCode, tree->zero, code);
        code.remove(code.size() - 1);

        // go to right
        code.add(1);
        // tree = tree->one;
        getCharCode(charAndCode, tree->one, code);
        code.remove(code.size() - 1);
    }
}

/**
 * Given a string and an encoding tree, encode the text using the tree
 * and return a Queue<Bit> of the encoded bit sequence.
 *
 * You can assume tree is a valid non-empty encoding tree and contains an
 * encoding for every character in the text.
 */
Queue<Bit> encodeText(EncodingTreeNode* tree, string text) {
    Queue<Bit> textCode;
    Map<char, Vector<Bit>> charsAndCodes;
    Vector<Bit> code;

    // get the map
    getCharCode(charsAndCodes, tree, code);
    for (char charcater : charsAndCodes) {
        cout << charcater << "-----" << charsAndCodes[charcater] << endl;
    }


    // according to the map, encode the text
    for (char character : text) {
        if (charsAndCodes.containsKey(character)) {
            Vector<Bit> charCode = charsAndCodes[character];
            for (Bit singleCode : charCode) {
                textCode.enqueue(singleCode);
            }
        }
    }

    return textCode;
}

/**
 * Flatten the given tree into a Queue<Bit> and Queue<char> in the manner
 * specified in the assignment writeup.
 *
 * You can assume the input Queues are empty on entry to this function.
 *
 * You can assume tree is a valid well-formed encoding tree.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 */
void flattenTree(EncodingTreeNode* tree, Queue<Bit>& treeBits, Queue<char>& treeLeaves) {
    /* TODO: Fill in the rest of this function. */
    if (tree == nullptr) {
        return;
    }

    if (isLeaf(tree)){ // is leaf, get the char and enqueue to treeleaves
        treeBits.enqueue(0);
        treeLeaves.enqueue(tree->ch);
        return;
    } else { // is connector, traversal down
        treeBits.enqueue(1);
        flattenTree(tree->zero, treeBits, treeLeaves);
        flattenTree(tree->one, treeBits, treeLeaves);
    }

}

/**
 * Compress the input text using Huffman coding, producing as output
 * an EncodedData containing the encoded message and encoding tree used.
 *
 * Reports an error if the message text does not contain at least
 * two distinct characters.
 *
 * TODO: Add any additional information to this comment that is necessary to describe
 * your implementation.
 */
EncodedData compress(string messageText) {
    EncodingTreeNode *root = buildHuffmanTree(messageText);
    Queue<Bit> messageCode = encodeText(root, messageText);

    Queue<Bit> treeBits;
    Queue<char> treeLeaves;
    flattenTree(root, treeBits, treeLeaves);

    EncodedData encodedMessageTree = {treeBits, treeLeaves, messageCode};

    deallocateTree(root);

    return encodedMessageTree;
}

//EncodingTreeNode* reference = createExampleTree(); // see diagram above
//EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
//EXPECT(areEqual(tree, reference));

//deallocateTree(reference);
//deallocateTree(tree);

/* * * * * * Testing Helper Functions Below This Point * * * * * */

EncodingTreeNode* createExampleTree() {
    /* Example encoding tree used in multiple test cases:
     *                *
     *              /   \
     *             T     *
     *                  / \
     *                 *   E
     *                / \
     *               R   S
     */
    /* TODO: Implement this helper function needed for testing. */
    EncodingTreeNode* r = new EncodingTreeNode('R');
    EncodingTreeNode* s = new EncodingTreeNode('S');
    EncodingTreeNode* e = new EncodingTreeNode('E');
    EncodingTreeNode* t = new EncodingTreeNode('T');

    EncodingTreeNode* rs = new EncodingTreeNode(r, s);
    EncodingTreeNode* rse = new EncodingTreeNode(rs, e);
    EncodingTreeNode* root = new EncodingTreeNode(t, rse);

    return root;
}

void deallocateTree(EncodingTreeNode* root) {
    /* TODO: Implement this helper function needed for testing. */
    if (root == nullptr) {
        return;
    }

    deallocateTree(root->zero);
    deallocateTree(root->one);
    delete root;
}


bool areEqual(EncodingTreeNode* a, EncodingTreeNode* b) {
    /* TODO: Implement this helper function needed for testing. */
    if (a == nullptr || b == nullptr) {
        return a == b;
    }

    return (a->ch == b->ch) && areEqual(a->zero, b->zero) && areEqual(a->one, b->one);

//    if (isLeaf(a) == false && isLeaf(b) == false) {
//        areEqual(a->zero, b->zero);
//        areEqual(a->one, b->one);
//    } else if (isLeaf(a) == true && isLeaf(b) == true) {
//        if (a->ch == b->ch) {
//            return true;
//        } else {
//            return false;
//        }
//    }
//    return false;
}

/* * * * * * Test Cases Below This Point * * * * * */

/* TODO: Write your own student tests. */









/* * * * * Provided Tests Below This Point * * * * */

PROVIDED_TEST("decodeText, small example encoding tree") {
    EncodingTreeNode* tree = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(decodeText(tree, messageBits), "E");

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(decodeText(tree, messageBits), "SET");

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1}; // STREETS
    EXPECT_EQUAL(decodeText(tree, messageBits), "STREETS");

    deallocateTree(tree);
}

PROVIDED_TEST("unflattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  treeBits   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeLeaves = { 'T', 'R', 'S', 'E' };
    EncodingTreeNode* tree = unflattenTree(treeBits, treeLeaves);

    EXPECT(areEqual(tree, reference));

    deallocateTree(tree);
    deallocateTree(reference);
}

PROVIDED_TEST("decompress, small example input") {
    EncodedData data = {
        { 1, 0, 1, 1, 0, 0, 0 }, // treeBits
        { 'T', 'R', 'S', 'E' },  // treeLeaves
        { 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1 } // messageBits
    };

    EXPECT_EQUAL(decompress(data), "TRESS");
}

PROVIDED_TEST("buildHuffmanTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    EncodingTreeNode* tree = buildHuffmanTree("STREETTEST");
    EXPECT(areEqual(tree, reference));

    deallocateTree(reference);
    deallocateTree(tree);
}

PROVIDED_TEST("encodeText, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above

    Queue<Bit> messageBits = { 1, 1 }; // E
    EXPECT_EQUAL(encodeText(reference, "E"), messageBits);

    messageBits = { 1, 0, 1, 1, 1, 0 }; // SET
    EXPECT_EQUAL(encodeText(reference, "SET"), messageBits);

    messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1 }; // STREETS
    EXPECT_EQUAL(encodeText(reference, "STREETS"), messageBits);

    deallocateTree(reference);
}

PROVIDED_TEST("flattenTree, small example encoding tree") {
    EncodingTreeNode* reference = createExampleTree(); // see diagram above
    Queue<Bit>  expectedBits   = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> expectedLeaves = { 'T', 'R', 'S', 'E' };

    Queue<Bit>  treeBits;
    Queue<char> treeLeaves;
    flattenTree(reference, treeBits, treeLeaves);

    EXPECT_EQUAL(treeBits,   expectedBits);
    EXPECT_EQUAL(treeLeaves, expectedLeaves);

    deallocateTree(reference);
}

PROVIDED_TEST("compress, small example input") {
    EncodedData data = compress("STREETTEST");
    Queue<Bit>  treeBits    = { 1, 0, 1, 1, 0, 0, 0 };
    Queue<char> treeChars   = { 'T', 'R', 'S', 'E' };
    Queue<Bit>  messageBits = { 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0 };

    EXPECT_EQUAL(data.treeBits, treeBits);
    EXPECT_EQUAL(data.treeLeaves, treeChars);
    EXPECT_EQUAL(data.messageBits, messageBits);
}

PROVIDED_TEST("Test end-to-end compress -> decompress") {
    Vector<string> inputs = {
        "HAPPY HIP HOP",
        "The job requires extra pluck and zeal from every young wage earner.",
        ":-) :-D XD <(^_^)>",
    };

    for (string input: inputs) {
        EncodedData data = compress(input);
        string output = decompress(data);

        EXPECT_EQUAL(output.size(), input.size());

        /* Don't clobber the output with a huge string if there's a mismatch. */
        EXPECT(input == output);
    }
}
