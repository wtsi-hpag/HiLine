/*
Copyright (c) 2020 Ed Harry, Wellcome Sanger Institute

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// inspired by https://noulin.net/wavlTree/file/wavlTree.c.html

struct
wavl_node
{
    wavl_node *parent;
    wavl_node *left;
    wavl_node *right;
    wavl_node *next;
    u64 value;
    u64 rank;
    u08 *tag;
    u64 tagLength;
};

struct
wavl_tree
{
    wavl_node *root;
    u64 numIntervals;
};

global_function
wavl_node *
WavlTreeNewNode(memory_arena *arena, wavl_node *parent, u64 value)
{
    wavl_node *node = PushStructP(arena, wavl_node);
    node->parent = parent;
    node->left = 0;
    node->right = 0;
    node->next = 0;
    node->value = value;
    node->rank = 0;

    return(node);
}

global_function
wavl_tree *
InitialiseWavlTree(memory_arena *arena)
{
    wavl_tree *tree = PushStructP(arena, wavl_tree);
    tree->numIntervals = 1;
    tree->root = WavlTreeNewNode(arena, 0, 0);
    
    return(tree);
}

global_function
u32
WavlTreeNeedToRotateLeftStrong(wavl_node *node)
{
    return((node->rank == 1 && !node->left) || (node->rank > node->left->rank + 1));
}

global_function
u32
WavlTreeNeedToRotateRightStrong(wavl_node *node)
{
    return((node->rank == 1 && !node->right) || (node->rank > node->right->rank + 1));
}

global_function
u32
WavlTreeNeedToRotateLeftWeak(wavl_node *node)
{
    return((!node->left) || (node->rank > node->left->rank + 1));
}

global_function
u32
WavlTreeNeedToRotateRightWeak(wavl_node *node)
{
    return((!node->right) || (node->rank > node->right->rank + 1));
}

global_function
void
WavlTreeRotateLeft(wavl_tree *tree, wavl_node *node)
{
    wavl_node *right = node->right;
    node->right = right->left;
    if (right->left)
    {
        right->left->parent = node;
    }
    right->parent = node->parent;
    if (!node->parent)
    {
        tree->root = right;
    }
    else if (node->parent->left == node)
    {
        node->parent->left = right;
    }
    else
    {
        node->parent->right = right;
    }
    right->left = node;
    node->parent = right;
}

global_function
void
WavlTreeRotateRight(wavl_tree *tree, wavl_node *node)
{
    wavl_node *left = node->left;
    node->left = left->right;
    if (left->right)
    {
        left->right->parent = node;
    }
    left->parent = node->parent;
    if (!node->parent)
    {
        tree->root = left;
    }
    else if (node->parent->right == node)
    {
        node->parent->right = left;
    }
    else
    {
        node->parent->left = left;
    }
    left->right = node;
    node->parent = left;
}

global_function
void
WavlTreeBalance(wavl_tree *tree, wavl_node *node)
{
    for (   wavl_node *parent = node->parent;
            parent && ((node->rank + 1) != parent->rank);
            node = parent, parent = node->parent, ++node->rank )
    {
        if (parent->left == node)
        {
            if (WavlTreeNeedToRotateRightStrong(parent))
            {
                if (WavlTreeNeedToRotateLeftWeak(node))
                {
                    --node->rank;
                    ++node->right->rank;
                    WavlTreeRotateLeft(tree, node);
                }
                --parent->rank;
                WavlTreeRotateRight(tree, parent);
                break;
            }
        }
        else if (WavlTreeNeedToRotateLeftStrong(parent))
        {
            if (WavlTreeNeedToRotateRightWeak(node))
            {
                --node->rank;
                ++node->left->rank;
                WavlTreeRotateRight(tree, node);
            }
            --parent->rank;
            WavlTreeRotateLeft(tree, parent);
            break;
        }
    }
}

global_function
void
WavlTreeInsertValue(memory_arena *arena, wavl_tree *tree, u64 value)
{
    wavl_node *node = tree->root;
    wavl_node *parent;
    do
    {
        if (value == node->value) return;
        parent = node;
        node = value < node->value ? node->left : node->right;
    } while (node);

    wavl_node *newNode = WavlTreeNewNode(arena, parent, value);
    ++tree->numIntervals;
    (value < parent->value ? parent->left : parent->right) = newNode;

    if (!parent->rank)
    {
        parent->rank = 1;
        WavlTreeBalance(tree, parent);
    }
}

global_function
wavl_node *
WavlTreeFindInterval(wavl_tree *tree, u64 value)
{
    wavl_node *result = 0;
    wavl_node *node = tree->root;
    do
    {
        if (value < node->value)
        {
            node = node->left;
        }
        else
        {
            result = node;
            node = node->right;
        }

    } while (node);

    return(result);
}

global_function
void
WavlTreeFreeze(wavl_node *node, wavl_node **prevNode)
{
    if (node) 
    { 
        WavlTreeFreeze(node->left, prevNode); 

        if (*prevNode) (*prevNode)->next = node;
        *prevNode = node;

        WavlTreeFreeze(node->right, prevNode); 
    } 
}

global_function
void
WavlTreeFreeze(wavl_tree *tree)
{
    wavl_node *tmp = 0;
    WavlTreeFreeze(tree->root, &tmp);
}

