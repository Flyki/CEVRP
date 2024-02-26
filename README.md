# CEVRP

Using BACO to solve the basic CEVRP. Please read "Jia, Y. H., Mei, Y., & Zhang, M. (2021). A Bilevel Ant Colony Optimization Algorithm for Capacitated Electric Vehicle Routing Problem. IEEE Transactions on Cybernetics." for the explanation of the algorithm.

Using CACO to solve the basic CEVRP with comparison of different encoding schemes. Please read "Jia, Y. H., Mei, Y., & Zhang, M. (2022). Confidence-based Ant Colony Optimization for Capacitated Electric Vehicle Routing Problem with Comparison of Different Encoding Schemes. IEEE Transactions on Evolutionary Computation." for the explanation of the algorithm.



> The main changes made by me is on three aspects:
>
> - Max-Evals stop criteria (aligns with the competition benchmark stop criteria)
> - Statistical function
> - Multithreading running

- Usage:

  1. First step - compile

     ```shell
     mkdir build
     cd build
     cmake ..
     make
     ```

  2. Second step - run

     ```shell
     ./build/CEVRP-Yahui cbaco E-n23-k3.evrp 1 1 1
     
     # Explanation
     # ./build/CEVRP-Yahui <algorithm: baco, cbaco> <problem_instance_filename> <pop_init> <confidence_based_selection> <representation>
     ```

     

  

