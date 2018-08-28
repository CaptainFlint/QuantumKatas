// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

namespace Quantum.Kata.Superposition
{
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;

    //////////////////////////////////////////////////////////////////
    // Welcome!
    //////////////////////////////////////////////////////////////////

    // "Superposition" quantum kata is a series of exercises designed 
    // to get you familiar with programming in Q#.
    // It covers the following topics:
    //  - basic single-qubit and multi-qubit gates,
    //  - superposition,
    //  - flow control and recursion in Q#.
    //
    // Each task is wrapped in one operation preceded by the description of the task.
    // Each task (except tasks in which you have to write a test) has a unit test associated with it,
    // which initially fails. Your goal is to fill in the blank (marked with // ... comment)
    // with some Q# code to make the failing test pass.
    //
    // The tasks are given in approximate order of increasing difficulty; harder ones are marked with asterisks. 

    // Task 1. Plus state
    // Input: a qubit in |0⟩ state (stored in an array of length 1).
    // Goal: create a |+⟩ state on this qubit (|+⟩ = (|0⟩ + |1⟩) / sqrt(2)).
    operation PlusState (qs : Qubit[]) : ()
    {
        body
        {
            // Hadamard gate H will convert |0⟩ state to |+⟩ state.
            // The first qubit of the array can be accessed as qs[0].
            // Type H(qs[0]);
            // Then rebuild the project and rerun the tests - T01_PlusState_Test should now pass!

            H(qs[0]);
        }
    }

    // Task 2. Minus state
    // Input: a qubit in |0⟩ state (stored in an array of length 1).
    // Goal: create a |-⟩ state on this qubit (|-⟩ = (|0⟩ - |1⟩) / sqrt(2)).
    operation MinusState (qs : Qubit[]) : ()
    {
        body
        {
            // In this task, as well as in all subsequent ones, you have to come up with the solution yourself.
            
            X(qs[0]);
            H(qs[0]);
        }
    }

    // Task 3*. Unequal superposition
    // Inputs:
    //      1) a qubit in |0⟩ state (stored in an array of length 1).
    //      2) angle alpha, in radians, represented as Double
    // Goal: create a cos(alpha) * |0⟩ + sin(alpha) * |1⟩ state on this qubit.
    operation UnequalSuperposition (qs : Qubit[], alpha : Double) : ()
    {
        body
        {
            // Hint: Experiment with rotation gates from Microsoft.Quantum.Primitive namespace.
            // Note that all rotation operators rotate the state by _half_ of its angle argument.

            Ry(alpha * 2.0, qs[0]);
        }
    }

    // Task 4. Bell state
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal: create a Bell state |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2) on these qubits.
    operation BellState (qs : Qubit[]) : ()
    {
        body
        {
            H(qs[0]);
            CNOT(qs[0], qs[1]);
        }
    }

    // Task 5. All Bell states
    // Inputs:
    //      1) two qubits in |00⟩ state (stored in an array of length 2)
    //      2) an integer index
    // Goal: create one of the Bell states based on the value of index:
    //       0: |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2)
    //       1: |Φ⁻⟩ = (|00⟩ - |11⟩) / sqrt(2)
    //       2: |Ψ⁺⟩ = (|01⟩ + |10⟩) / sqrt(2)
    //       3: |Ψ⁻⟩ = (|01⟩ - |10⟩) / sqrt(2)
    operation AllBellStates (qs : Qubit[], index : Int) : ()
    {
        body
        {
            if (index == 0) {
                H(qs[0]);
                CNOT(qs[0], qs[1]);
            }
            elif (index == 1) {
                X(qs[0]);
                H(qs[0]);
                CNOT(qs[0], qs[1]);
            }
            elif (index == 2) {
                H(qs[0]);
                CNOT(qs[0], qs[1]);
                X(qs[1]);
            }
            elif (index == 3) {
                X(qs[0]);
                H(qs[0]);
                CNOT(qs[0], qs[1]);
                X(qs[1]);
            }
        }
    }

    // Task 6. Greenberger–Horne–Zeilinger state
    // Input: N qubits in |0...0⟩ state.
    // Goal: create a GHZ state (|0...0⟩ + |1...1⟩) / sqrt(2) on these qubits.
    operation GHZ_State (qs : Qubit[]) : ()
    {
        body
        {
            // Hint: N can be found as Length(qs).

            H(qs[0]);
            for (index in 1 .. Length(qs) - 1) {
                CNOT(qs[0], qs[index]);
            }
        }
    }

    // Task 7. Superposition of all basis vectors
    // Input: N qubits in |0...0⟩ state.
    // Goal: create an equal superposition of all basis vectors from |0...0⟩ to |1...1⟩
    // (i.e. state (|0...0⟩ + ... + |1...1⟩) / sqrt(2^N) ).
    operation AllBasisVectorsSuperposition (qs : Qubit[]) : ()
    {
        body
        {
            for (index in 0 .. Length(qs) - 1) {
                H(qs[index]);
            }
        }
    }

    // Task 8. Superposition of |0...0⟩ and given bit string
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) bit string represented as Bool[]
    // Goal: create an equal superposition of |0...0⟩ and basis state given by the bit string.
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // You are guaranteed that the qubit array and the bit string have the same length.
    // You are guaranteed that the first bit of the bit string is true.
    // Example: for bit string = [true; false] the qubit state required is (|00⟩ + |10⟩) / sqrt(2).
    operation ZeroAndBitstringSuperposition (qs : Qubit[], bits : Bool[]) : ()
    {
        body
        {
            // The following lines enforce the constraints on the input that you are given.
            // You don't need to modify them. Feel free to remove them, this won't cause your code to fail.
            AssertIntEqual(Length(bits), Length(qs), "Arrays should have the same length");
            AssertBoolEqual(bits[0], true, "First bit of the input bit string should be set to true");

            H(qs[0]);
            for (index in 1 .. Length(qs) - 1) {
                if (bits[index]) {
                    CNOT(qs[0], qs[index]);
                }
            }
        }
    }

    // Task 9. Superposition of two bit strings
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) two bit string represented as Bool[]s
    // Goal: create an equal superposition of two basis states given by the bit strings.
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // Example: for bit strings [false; true; false] and [false; false; true] 
    // the qubit state required is (|010⟩ + |001⟩) / sqrt(2).
    // You are guaranteed that the both bit strings have the same length as the qubit array,
    // and that the bit strings will differ in at least one bit.
    operation TwoBitstringSuperposition (qs : Qubit[], bits1 : Bool[], bits2 : Bool[]) : ()
    {
        body
        {
            mutable index1 = 0;
            // Process to the first difference
            repeat {
                if (bits1[index1] && bits2[index1]) {
                    X(qs[index1]);
                }
            } until (bits1[index1] != bits2[index1])
            fixup {
                set index1 = index1 + 1;
            }
            let ctrl = index1;
            let ctrl_one = bits1[ctrl]; // Controller bit has true in the first array
            H(qs[ctrl]);
            // Now we have superposition, proceed accordingly
            for (index2 in (index1 + 1) .. (Length(qs) - 1)) {
                if (bits1[index2] == bits2[index2]) {
                    if (bits1[index2]) {
                        X(qs[index2]);
                    }
                }
                else {
                    CNOT(qs[ctrl], qs[index2]);
                    if (bits1[index2] != ctrl_one) {
                        X(qs[index2]);
                    }
                }
            }
        }
    }

    // Task 10**. W state on 2^k qubits
    // Input: N = 2^k qubits in |0...0⟩ state.
    // Goal: create a W state (https://en.wikipedia.org/wiki/W_state) on these qubits.
    // W state is an equal superposition of all basis states on N qubits of Hamming weight 1.
    // Example: for N = 4, W state is (|1000⟩ + |0100⟩ + |0010⟩ + |0001⟩) / 2.
    operation WState_PowerOfTwo (qs : Qubit[]) : ()
    {
        body
        {
            // Hint: you can use Controlled functor to perform arbitrary controlled gates.

            let n = Length(qs);
            mutable k = 0;
            mutable N = 1;
            repeat {
            } until (N >= n)
            fixup {
                set k = k + 1;
                set N = N * 2;
            }
            if (k == 0) {
                // N = 1, special case: just return |1⟩
                X(qs[0]);
            }
            else {
                // Allocate additional k qubits that will act as generators for the final states
                using (qs_aux = Qubit[k]) {
                    // First, prepare an equal superposition of all possible states of qs_aux (N = 2^k states)
                    for (j in 0 .. k - 1) {
                        H(qs_aux[j]);
                    }
                    // Now for each of these states let's set into |1⟩ a unique input qubit.
                    // Matching rule: the ordinal index of the qs_aux state corresponds to its binary representation
                    // (e.g. with k = 4, the state No.6 will be |0110⟩). We will set to |1⟩ the input qubit
                    // whose index is equal to that state's index.
                    // The setting itself is performed by transforming the each aux state in turn into |1...1⟩
                    // by NOTting all the 0-qubits, calling controlled not from all the aux qubits onto the target
                    // qubit, then reverting the aux state.
                    //
                    // Example: N = 4, k = 2
                    // 1. Allocate 2 aux qubits |00⟩, full state: |0000.00⟩
                    // 2. Make superposition: |0000.00⟩ + |0000.01⟩ + |0000.10⟩ + |0000.11⟩
                    // 3. Set target qubits:
                    //   a) index = 0, the corresponding aux state is |00⟩, so to make it |11⟩ we apply NOT to both aux qubits:
                    //    |0000.00⟩ + |0000.01⟩ + |0000.10⟩ + |0000.11⟩ -> |0000.11⟩ + |0000.10⟩ + |0000.01⟩ + |0000.00⟩
                    //   b) Switch 0 qubit where aux state is |11⟩ using CCNOT:
                    //    |1000.11⟩ + |0000.10⟩ + |0000.01⟩ + |0000.00⟩
                    //   c) Revert aux state by repeating NOTs from 3.a:
                    //    |1000.00⟩ + |0000.01⟩ + |0000.10⟩ + |0000.11⟩
                    //   d) Repeat steps a-c for index in (1..3):
                    //    |1000.00⟩ + |0100.01⟩ + |0010.10⟩ + |0001.11⟩
                    // 4. After that reset the aux qubits separately in each state, using CNOT with the input qubits as controllers:
                    //    |1000.00⟩ + |0100.00⟩ + |0010.00⟩ + |0001.00⟩
                    //
                    // OK, let's rock and roll.
                    mutable bv = 0;
                    for (i in 0 .. n - 1) {     // i is the index of the input qubit / aux state
                        // 1. Switch the next aux state into |1...1⟩
                        set bv = 1;
                        for (j in 0 .. k - 1) { // j is the index of the aux qubit
                            if ((i &&& bv) == 0) {
                                X(qs_aux[j]);
                            }
                            set bv = bv * 2;
                        }
                        // 2. Set the corresponding input qubit into |1⟩
                        (Controlled(X))(qs_aux, (qs[i]));
                        // 3. Restore the aux state
                        set bv = 1;
                        for (j in 0 .. k - 1) {
                            if ((i &&& bv) == 0) {
                                X(qs_aux[j]);
                            }
                            set bv = bv * 2;
                        }
                    }
                    // Now we need to reset all the aux qubits into zero state without affecting the
                    // main input qubits that are now entangled with them.
                    // Again, just treat the aux states as binary representations, but this time reset
                    // the 1-bits using CNOT from the main qubit (it's guaranteed to be the only one in |1⟩)
                    for (i in 0 .. n - 1) {
                        set bv = 1;
                        for (j in 0 .. k - 1) {
                            if ((i &&& bv) != 0) {
                                CNOT(qs[i], qs_aux[j]);
                            }
                            set bv = bv * 2;
                        }
                    }
                }
            }
        }
    }

    // Task 11**. W state on arbitrary number of qubits
    // Input: N qubits in |0...0⟩ state (N is not necessarily a power of 2).
    // Goal: create a W state (https://en.wikipedia.org/wiki/W_state) on these qubits.
    // W state is an equal superposition of all basis states on N qubits of Hamming weight 1.
    // Example: for N = 3, W state is (|100⟩ + |010⟩ + |001⟩) / sqrt(3).
    operation WState_Arbitrary (qs : Qubit[]) : ()
    {
        body
        {
            let n = Length(qs);
            mutable k = 0;
            mutable N = 1;
            repeat {
            } until (N >= n)
            fixup {
                set k = k + 1;
                set N = N * 2;
            }
            Message($"len={n}, N={N}, k={k}");
            if (N == n) {
                // For power of two call the already implemented function
                WState_PowerOfTwo(qs);
            }
            else {
                // Supplemental qubits to make a total of N = 2^k qubits
                using (qs_sup = Qubit[N - n]) {
                    // Prepare the W-state for N qubits.
                    // Examples for n = 5 and n = 6 (N = 8):
                    //   |10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000> + |00000.100> + |00000.010> + |00000.001>
                    //   |100000.00> + |010000.00> + |001000.00> + |000100.00> + |000010.00> + |000001.00> + |000000.10> + |000000.01>
                    WState_PowerOfTwo(qs + qs_sup);

                    // One-by-one, merge the supplemental qubits into the first state:
                    //   a) CNOT for setting the first qubit in the supplemental state;
                    //   b) controlled Ry for merging the first state with the supplemental one, so that supplemental qubit was |0>,
                    //      the rotation angle equals -2 * arctg(1/√m) where √m is the amplitude of the |0> state.
                    // Chain for n = 5:
                    //   |10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000> + |10000.100> + |00000.010> + |00000.001>
                    //   √2|10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000> + |00000.010> + |00000.001>
                    //   √2|10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000> + |10000.010> + |00000.001>
                    //   √3|10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000> + |00000.001>
                    //   √3|10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000> + |10000.001>
                    //   √4|10000.000> + |01000.000> + |00100.000> + |00010.000> + |00001.000>
                    // Final result for n = 6:
                    //   √3|100000.00> + |010000.00> + |001000.00> + |000100.00> + |000010.00> + |000001.00>
                    for (i in 0 .. N - n - 1) {
                        CNOT(qs_sup[i], qs[0]);
                        (Controlled(Ry))([qs[0]], (-2.0 * ArcTan(1.0 / Sqrt(ToDouble(i + 1))), qs_sup[i]));
                    }
                    // Now supplemental qubits are reset, we can get rid of them.
                }
                // Current state:
                //   √4|10000> + |01000> + |00100> + |00010> + |00001>
                //   √3|100000> + |010000> + |001000> + |000100> + |000010> + |000001>
                // Using Ry...
                for (i in 0 .. N - n - 1) {
                    //(controlled(Ry))([qs[i]], qs_sup[i]);
                }
            }
        }
    }
}
