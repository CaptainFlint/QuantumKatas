﻿// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

namespace Quantum.Kata.DeutschJozsaAlgorithm
{
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;

    //////////////////////////////////////////////////////////////////
    // Welcome!
    //////////////////////////////////////////////////////////////////

    // "Deutsch-Jozsa algorithm" quantum kata is a series of exercises designed 
    // to get you familiar with programming in Q#.
    // It covers the following topics:
    //  - writing oracles (quantum operations which implement certain classical functions),
    //  - Bernstein-Vazirani algorithm for recovering the parameters of a scalar product function,
    //  - Deutsch-Jozsa algorithm for recognizing a function as constant or balanced, and
    //  - writing tests in Q#.
    //
    // Each task is wrapped in one operation preceded by the description of the task.
    // Each task (except tasks in which you have to write a test) has a unit test associated with it,
    // which initially fails. Your goal is to fill in the blank (marked with // ... comment)
    // with some Q# code to make the failing test pass.

    //////////////////////////////////////////////////////////////////
    // Part I. Oracles
    //////////////////////////////////////////////////////////////////

    // In this section you will implement oracles defined by classical functions using the following rules:
    //  - a function f(𝑥₀, …, 𝑥ₙ₋₁) with N bits of input x = (𝑥₀, …, 𝑥ₙ₋₁) and 1 bit of output y
    //    defines an oracle which acts on N input qubits and 1 output qubit.
    //  - the oracle effect on qubits in computational basis states is defined as follows:
    //    |x⟩ |y⟩ -> |x⟩ |y ⊕ f(x)⟩   (⊕ is addition modulo 2)
    //  - the oracle effect on qubits in superposition is defined following the linearity of quantum operations.
    //  - the oracle must act properly on qubits in all possible input states.

    // Task 1.1. f(x) = 0
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    operation Oracle_Zero (x : Qubit[], y : Qubit) : ()
    {
        body
        {
            // Since f(x) = 0 for all values of x, |y ⊕ f(x)⟩ = |y⟩.
            // This means that the operation doesn't need to do any transformation to the inputs.
            // Build the project and run the tests to see that T01_Oracle_Zero_Test test passes.
        }
    }

    // Task 1.2. f(x) = 1
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    operation Oracle_One (x : Qubit[], y : Qubit) : ()
    {
        body
        {
            // Since f(x) = 1 for all values of x, |y ⊕ f(x)⟩ = |y ⊕ 1⟩ = |NOT y⟩.
            // This means that the operation needs to flip qubit y (i.e. transform |0⟩ to |1⟩ and vice versa).

            X(y);
        }
    }

    // Task 1.3. f(x) = xₖ (the value of k-th qubit)
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    //      3) 0-based index of the qubit from input register (0 <= k < N)
    // Goal: transform state |x, y⟩ into state |x, y ⊕ xₖ⟩ (⊕ is addition modulo 2).
    operation Oracle_Kth_Qubit (x : Qubit[], y : Qubit, k : Int) : ()
    {
        body
        {
            // The following line enforces the constraints on the value of k that you are given.
            // You don't need to modify it. Feel free to remove it, this won't cause your code to fail.
            AssertBoolEqual(0 <= k && k < Length(x), true, "k should be between 0 and N-1, inclusive");

            CNOT(x[k], y);
        }
    }

    // Task 1.4. f(x) = 1 if x has odd number of 1s, and 0 otherwise
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    operation Oracle_OddNumberOfOnes (x : Qubit[], y : Qubit) : ()
    {
        body
        {
            // Hint: f(x) can be represented as x_0 ⊕ x_1 ⊕ ... ⊕ x_(N-1)

            for (i in 0 .. Length(x) - 1) {
                CNOT(x[i], y);
            }
        }
    }

    // Task 1.5. f(x) = Σᵢ 𝑟ᵢ 𝑥ᵢ modulo 2 for a given bit vector r (scalar product function)
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    //      3) a bit vector of length N represented as Int[]
    // You are guaranteed that the qubit array and the bit vector have the same length.
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    //
    // Note: the functions featured in tasks 1.1, 1.3 and 1.4 are special cases of this function.
    operation Oracle_ProductFunction (x : Qubit[], y : Qubit, r : Int[]) : ()
    {
        body
        {
            // The following line enforces the constraint on the input arrays.
            // You don't need to modify it. Feel free to remove it, this won't cause your code to fail.
            AssertIntEqual(Length(x), Length(r), "Arrays should have the same length");

            for (i in 0 .. Length(x) - 1) {
                if (r[i] != 0) {
                    CNOT(x[i], y);
                }
            }
        }
    }

    // Task 1.6. f(x) = Σᵢ (𝑟ᵢ 𝑥ᵢ + (1 - 𝑟ᵢ)(1 - 𝑥ᵢ)) modulo 2 for a given bit vector r
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    //      3) a bit vector of length N represented as Int[]
    // You are guaranteed that the qubit array and the bit vector have the same length.
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    operation Oracle_ProductWithNegationFunction (x : Qubit[], y : Qubit, r : Int[]) : ()
    {
        body
        {
            // The following line enforces the constraint on the input arrays.
            // You don't need to modify it. Feel free to remove it, this won't cause your code to fail.
            AssertIntEqual(Length(x), Length(r), "Arrays should have the same length");

            for (i in 0 .. Length(x) - 1) {
                if (r[i] != 0) {
                    CNOT(x[i], y);
                }
                else {
                    X(x[i]);
                    CNOT(x[i], y);
                    X(x[i]);
                }
            }
        }
    }

    // Task 1.7. f(x) = Σᵢ 𝑥ᵢ + (1 if prefix of x is equal to the given bit vector, and 0 otherwise) modulo 2
    // Inputs: 
    //      1) N qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    //      3) a bit vector of length P represented as Int[] (1 <= P <= N)
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    // 
    // A prefix of length k of a state |x⟩ = |x₁, ..., xₙ⟩ is the state of its first k qubits |x₁, ..., xₖ⟩.
    // For example, a prefix of length 2 of a state |0110⟩ is 01.
    operation Oracle_HammingWithPrefix (x : Qubit[], y : Qubit, prefix : Int[]) : ()
    {
        body
        {
            // The following line enforces the constraint on the input arrays.
            // You don't need to modify it. Feel free to remove it, this won't cause your code to fail.
            let P = Length(prefix);
            AssertBoolEqual(1 <= P && P <= Length(x), true, "P should be between 1 and N, inclusive");

            // Hint: the first part of the function is the same as in task 1.4

            for (i in 0 .. Length(x) - 1) {
                CNOT(x[i], y);
            }

            // Hint: you can use Controlled functor to perform multicontrolled gates
            // (gates with multiple control qubits).

            for (i in 0 .. P - 1) {
                if (prefix[i] == 0) {
                    X(x[i]);
                }
            }
            (Controlled(X))(x[0 .. P - 1], (y));
            for (i in 0 .. P - 1) {
                if (prefix[i] == 0) {
                    X(x[i]);
                }
            }
        }
    }

    // Task 1.8*. f(x) = 1 if x has two or three bits (out of three) set to 1, and 0 otherwise  (majority function)
    // Inputs: 
    //      1) 3 qubits in arbitrary state |x⟩ (input register)
    //      2) a qubit in arbitrary state |y⟩ (output qubit)
    // Goal: transform state |x, y⟩ into state |x, y ⊕ f(x)⟩ (⊕ is addition modulo 2).
    operation Oracle_MajorityFunction (x : Qubit[], y : Qubit) : ()
    {
        body
        {
            // The following line enforces the constraint on the input array.
            // You don't need to modify it. Feel free to remove it, this won't cause your code to fail.
            AssertBoolEqual(3 == Length(x), true, "x should have exactly 3 qubits");
            
            // Hint: represent f(x) in terms of AND and ⊕ operations

            // x = |111⟩
            (Controlled(X))(x, (y));
            // x = |011⟩
            X(x[0]);
            (Controlled(X))(x, (y));
            X(x[0]);
            // x = |101⟩
            X(x[1]);
            (Controlled(X))(x, (y));
            X(x[1]);
            // x = |110⟩
            X(x[2]);
            (Controlled(X))(x, (y));
            X(x[2]);
        }
    }


    //////////////////////////////////////////////////////////////////
    // Part II. Bernstein-Vazirani Algorithm
    //////////////////////////////////////////////////////////////////

    // Task 2.1. State preparation for Bernstein-Vazirani algorithm
    // Inputs:
    //      1) N qubits in |0⟩ state (query register)
    //      2) a qubit in |0⟩ state (answer register)
    // Goal:
    //      1) create an equal superposition of all basis vectors from |0...0⟩ to |1...1⟩ on query register
    //         (i.e. state (|0...0⟩ + ... + |1...1⟩) / sqrt(2^N) )
    //      2) create |-⟩ state (|-⟩ = (|0⟩ - |1⟩) / sqrt(2)) on answer register
    operation BV_StatePrep (query : Qubit[], answer : Qubit) : ()
    {
        body
        {
            for (i in 0 .. Length(query) - 1) {
                H(query[i]);
            }
            X(answer);
            H(answer);
        }
        adjoint auto;
    }

    // Task 2.2. Bernstein-Vazirani algorithm implementation
    // Inputs:
    //      1) the number of qubits in the input register N for the function f
    //      2) a quantum operation which implements the oracle |x⟩|y⟩ -> |x⟩|y ⊕ f(x)⟩, where
    //         x is N-qubit input register, y is 1-qubit answer register, and f is a Boolean function
    // You are guaranteed that the function f implemented by the oracle is a scalar product function
    // (can be represented as f(𝑥₀, …, 𝑥ₙ₋₁) = Σᵢ 𝑟ᵢ 𝑥ᵢ modulo 2 for some bit vector r = (𝑟₀, …, 𝑟ₙ₋₁)).
    // You have implemented the oracle implementing the scalar product function in task 1.5.
    // Output:
    //      A bit vector r reconstructed from the function
    //
    // Note: a trivial approach is to call the oracle N times: 
    //       |10...0⟩|0⟩ = |10...0⟩|r₀⟩, |010...0⟩|0⟩ = |010...0⟩|r₁⟩ and so on.
    // Quantum computing allows to perform this task in just one call to the oracle; try to implement this algorithm.
    operation BV_Algorithm (N : Int, Uf : ((Qubit[], Qubit) => ())) : Int[]
    {
        body
        {
            // Declare a Bool array in which the result will be stored;
            // the array has to be mutable to allow updating its elements.
            mutable r = new Int[N];

            using (xy = Qubit[N + 1]) {
                // Use the preparation from the task 2.1
                BV_StatePrep(xy[0 .. N - 1], xy[N]);
                Uf(xy[0 .. N - 1], xy[N]);
                // Now we have a state similar to what we needed to distinguish in Measurements, Task 1.10
                for (i in 0 .. N) {
                    H(xy[i]);
                }
                for (i in 0 .. N - 1) {
                    if (M(xy[i]) == One) {
                        set r[i] = 1;
                    }
                    else {
                        set r[i] = 0;
                    }
                }
                ResetAll(xy);
            }

            return r;
        }
    }

    // Auxiliary function to run tests on operations with bit array as input
    operation ArrayOperator_TestCaller (TestProc : ((Int[]) => ())) : ()
    {
        body
        {
            let MaxN_Full = 8;          // For values up to this we enumerate all possible combinations
            let MaxN_Random = 16;       // Then, up to this one just do some randoms
            let RandomIterations = 100; // Amount of random tests
            for (N in 1 .. MaxN_Full) {
                // Iterations with full set of combinations
                for (iter in 0 .. 2 ^ N - 1) {
                    mutable r = new Int[N];
                    // Fill in the `r' with binary representation of the value `iter'
                    mutable v = iter;
                    for (i in 0 .. N - 1) {
                        set r[i] = v % 2;
                        set v = v / 2;
                    }
                    TestProc(r);
                }
            }
            for (N in MaxN_Full + 1 .. MaxN_Random) {
                // Iterations with random arrays
                for (iter in 1 .. RandomIterations) {
                    mutable r = new Int[N];
                    for (i in 0 .. N - 1) {
                        set r[i] = RandomInt(2);
                    }
                    TestProc(r);
                }
            }
        }
    }

    // Basic check for BV_Test for an arbitrary array
    operation BV_TestAssert (r : Int[]) : ()
    {
        body
        {
            AssertIntArrayEqual(r, BV_Algorithm(Length(r), Oracle_ProductFunction(_, _, r)),
                        "Inconsistency between scalar product oracle and BV algorithm result");
        }
    }

    // Task 2.3. Testing Bernstein-Vazirani algorithm
    // Goal: use your implementation of Bernstein-Vazirani algorithm from task 2.2 to figure out 
    // what bit vector the scalar product function oracle from task 1.5 was using.
    // As a reminder, this oracle creates an operation f(x) = Σᵢ 𝑟ᵢ 𝑥ᵢ modulo 2 for a given bit vector r,
    // and Bernstein-Vazirani algorithm recovers that bit vector given the operation.
    operation BV_Test () : ()
    {
        body
        {
            // Hint: use Oracle_ProductFunction to implement the scalar product function oracle passed to BV_Algorithm.
            // Since Oracle_ProductFunction takes three arguments (Qubit[], Qubit and Int[]), 
            // and the operation passed to BV_Algorithm must take two arguments (Qubit[] and Qubit), 
            // you need to use partial application to fix the third argument (a specific value of a bit vector). 
            // 
            // You might want to use something like the following:
            // let oracle = Oracle_ProductFunction(_, _, [...your bit vector here...]);

            // Hint: use AssertIntArrayEqual function to assert that the return value of BV_Algorithm operation 
            // matches the expected value (i.e. the bit vector passed to Oracle_ProductFunction).

            // BV_Test appears in the list of unit tests for the solution; run it to verify your code.

            ArrayOperator_TestCaller(BV_TestAssert);
        }
    }


    //////////////////////////////////////////////////////////////////
    // Part III. Deutsch-Jozsa Algorithm
    //////////////////////////////////////////////////////////////////

    // Task 3.1. Deutsch-Jozsa algorithm implementation
    // Inputs:
    //      1) the number of qubits in the input register N for the function f
    //      2) a quantum operation which implements the oracle |x⟩|y⟩ -> |x⟩|y ⊕ f(x)⟩, where
    //         x is N-qubit input register, y is 1-qubit answer register, and f is a Boolean function
    // You are guaranteed that the function f implemented by the oracle is either 
    // constant (returns 0 on all inputs or 1 on all inputs) or 
    // balanced (returns 0 on exactly one half of the input domain and 1 on the other half).
    // Output:
    //      true if the function f is constant 
    //      false if the function f is balanced
    //
    // Note: a trivial approach is to call the oracle multiple times: 
    //       if the values for more than half of the possible inputs are the same, the function is constant.
    // Quantum computing allows to perform this task in just one call to the oracle; try to implement this algorithm.
    operation DJ_Algorithm (N : Int, Uf : ((Qubit[], Qubit) => ())) : Bool
    {
        body
        {
            // Declare Bool variable in which the result will be accumulated;
            // this variable has to be mutable to allow updating it.
            mutable isConstantFunction = true;

            // Hint: even though Deutsch-Jozsa algorithm operates on a wider class of functions
            // than Bernstein-Vazirani (i.e. functions which can not be represented as a scalar product, such as f(x) = 1), 
            // it can be expressed as running Bernstein-Vazirani algorithm
            // and then post-processing the return value classically

            using (xy = Qubit[N + 1]) {
                // Use the preparation from the task 2.1
                BV_StatePrep(xy[0 .. N - 1], xy[N]);
                Uf(xy[0 .. N - 1], xy[N]);
                for (i in 0 .. N) {
                    H(xy[i]);
                }
                // **I did not prove this, really, but there are very good indications that the following is correct.**
                // Constant function turns this into |0...0⟩|1⟩ for 0, or -|0...0⟩|1⟩ for 1, which is really the same.
                // Other functions will have at least one 1 in the x.
                for (i in 0 .. N - 1) {
                    if (M(xy[i]) == One) {
                        set isConstantFunction = false;
                    }
                }
                ResetAll(xy);
            }

            return isConstantFunction;
        }
    }

    // Basic check for DJ_Test with Oracle_ProductFunction for an arbitrary array
    operation DJ15_TestAssert (r : Int[]) : ()
    {
        body
        {
            let N = Length(r);
            mutable expected = true;
            for (i in 0 .. N - 1) {
                if (r[i] != 0) {
                    set expected = false;
                }
            }
            AssertBoolEqual(DJ_Algorithm(N, Oracle_ProductFunction(_, _, r)), expected,
                "Oracle_ProductFunction was detected erroneously");
        }
    }

    // Basic check for DJ_Test with Oracle_ProductFunction for an arbitrary array
    operation DJ16_TestAssert (r : Int[]) : ()
    {
        body
        {
            let N = Length(r);
            AssertBoolEqual(DJ_Algorithm(N, Oracle_ProductWithNegationFunction(_, _, r)), false,
                "Oracle_ProductWithNegationFunction was detected erroneously");
        }
    }

    // Task 3.2. Testing Deutsch-Jozsa algorithm
    // Goal: use your implementation of Deutsch-Jozsa algorithm from task 3.1 to test 
    // each of the oracles you've implemented in part I for being constant or balanced.
    operation DJ_Test () : ()
    {
        body
        {
            // Hint: you will need to use partial application to test Oracle_Kth_Qubit and Oracle_ParityFunction;
            // see task 2.3 for a description of how to do that.

            // Hint: use AssertBoolEqual function to assert that the return value of DJ_Algorithm operation matches the expected value

            // DJ_Test appears in the list of unit tests for the solution; run it to verify your code.

            // 1.1 Oracle_Zero
            AssertBoolEqual(DJ_Algorithm(1, Oracle_Zero), true, "Oracle_Zero was erroneously determined as balanced");

            // 1.2 Oracle_One
            AssertBoolEqual(DJ_Algorithm(1, Oracle_One), true, "Oracle_One was erroneously determined as balanced");

            // 1.3 Oracle_Kth_Qubit
            // Trying lengths of 2 to 16
            for (N in 2 .. 16) {
                for (k in 0 .. N - 1) {
                    let oracle = Oracle_Kth_Qubit(_, _, k);
                    AssertBoolEqual(DJ_Algorithm(N, oracle), false, "Oracle_Kth_Qubit was erroneously determined as constant");
                }
            }

            // 1.4 Oracle_OddNumberOfOnes
            AssertBoolEqual(DJ_Algorithm(1, Oracle_OddNumberOfOnes), false, "Oracle_OddNumberOfOnes was erroneously determined as constant");

            // 1.5 Oracle_ProductFunction
            ArrayOperator_TestCaller(DJ15_TestAssert);

            // 1.6 Oracle_ProductWithNegationFunction
            ArrayOperator_TestCaller(DJ16_TestAssert);

            // 1.7 Oracle_HammingWithPrefix
            // Separately check N == 1, P == 1, function is constant
            AssertBoolEqual(DJ_Algorithm(1, Oracle_HammingWithPrefix(_, _, [0])), true,
                $"Oracle_HammingWithPrefix was erroneously determined as balanced for N == 1");
            AssertBoolEqual(DJ_Algorithm(1, Oracle_HammingWithPrefix(_, _, [1])), true,
                $"Oracle_HammingWithPrefix was erroneously determined as balanced for N == 1");
            // For other combinations, if P < N then the function is balanced;
            // if P == N then it's neither constant, nor balanced, so exluding this case.
            let MaxN = 16;              // Maximum qubit array length to test
            let MaxP_Full = 7;          // For values up to this we enumerate all possible prefix combinations
            let RandomIterations = 100; // Amount of random tests for the longer prefixes
            for (N in 2 .. MaxN) {
                // Maximum possible value of P
                mutable MaxP = N - 1;
                mutable P_Full = MaxP_Full;
                if (P_Full > MaxP) {
                    set P_Full = MaxP;
                }
                // Iterations with full set of combinations
                for (P in 1 .. P_Full) {
                    for (iter in 0 .. 2 ^ P - 1) {
                        mutable r = new Int[P];
                        // Fill in the `r' with binary representation of the value `iter'
                        mutable v = iter;
                        for (i in 0 .. P - 1) {
                            set r[i] = v % 2;
                            set v = v / 2;
                        }
                        AssertBoolEqual(DJ_Algorithm(N, Oracle_HammingWithPrefix(_, _, r)), N == 1,
                            $"Oracle_HammingWithPrefix was erroneously determined as constant, params: N={N}, P={P}, r={r}");
                    }
                }
                // Iterations with random arrays
                for (P in P_Full + 1 .. MaxP) {
                    for (iter in 1 .. RandomIterations) {
                        mutable r = new Int[P];
                        for (i in 0 .. P - 1) {
                            set r[i] = RandomInt(2);
                        }
                        AssertBoolEqual(DJ_Algorithm(N, Oracle_HammingWithPrefix(_, _, r)), false,
                            $"Oracle_HammingWithPrefix was erroneously determined as constant, params: N={N}, P={P}, r={r}");
                    }
                }
            }

            // 1.8 Oracle_MajorityFunction
            AssertBoolEqual(DJ_Algorithm(3, Oracle_MajorityFunction), false, "Oracle_MajorityFunction was erroneously determined as constant");
        }
    }


    //////////////////////////////////////////////////////////////////
    // Part IV. Come up with your own algorithm!
    //////////////////////////////////////////////////////////////////

    // Task 4.1. Reconstruct the oracle from task 1.6
    // Inputs:
    //      1) the number of qubits in the input register N for the function f
    //      2) a quantum operation which implements the oracle |x⟩|y⟩ -> |x⟩|y ⊕ f(x)⟩, where
    //         x is N-qubit input register, y is 1-qubit answer register, and f is a Boolean function
    // You are guaranteed that the function f implemented by the oracle can be represented as
    // f(𝑥₀, …, 𝑥ₙ₋₁) = Σᵢ (𝑟ᵢ 𝑥ᵢ + (1 - 𝑟ᵢ)(1 - 𝑥ᵢ)) modulo 2 for some bit vector r = (𝑟₀, …, 𝑟ₙ₋₁).
    // You have implemented the oracle implementing this function in task 1.6.
    // Output:
    //      A bit vector r which generates the same oracle as the one you are given
    operation Noname_Algorithm (N : Int, Uf : ((Qubit[], Qubit) => ())) : Int[]
    {
        body
        {
            // Declare a Bool array in which the result will be stored;
            // the array has to be mutable to allow updating its elements.
            mutable r = new Int[N];

            // ...

            return r;
        }
    }
}
