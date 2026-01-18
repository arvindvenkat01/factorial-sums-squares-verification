# Arvind Venkat | MIT License | Verified on Google Colab | Jan 18, 2026

import time
from gmpy2 import is_square, mpz


# ==========================================
# CONFIGURATION
# ==========================================
MAX_N = 40             
START_KERNELS = [ [] ]  # Start from empty set
# ==========================================


class UltraFastSearcher:
    def __init__(self, max_n):
        self.max_n = max_n
        self.results = []
        
        print(f"--- PRE-CALCULATING TABLES (N={max_n}) ---")
        
        # 1. Factorials
        self.factorials = []
        curr = mpz(1)
        for i in range(1, max_n + 2):
            curr *= i
            if i <= max_n:
                self.factorials.append(curr) 


        # 2. Primes Setup
        all_primes = self.get_primes(500)
        
        self.tree_primes = []
        self.leaf_primes = []
        
        for p in all_primes:
            if p <= max_n + 2:
                self.tree_primes.append(p)
            else:
                self.leaf_primes.append(p)
                
        # 3. QR Lookups
        self.qr_lookups = {}
        for p in all_primes:
            is_qr = [False] * p
            for x in range(p):
                is_qr[pow(x, 2, p)] = True
            self.qr_lookups[p] = is_qr
            
        # 4. Modulo Lookup Tables
        self.fact_mods = []
        for idx in range(max_n):
            mods = {}
            val = self.factorials[idx]
            for p in self.tree_primes:
                mods[p] = int(val % p)
            self.fact_mods.append(mods)


        # 5. Check Map
        self.check_map = [[] for _ in range(max_n + 1)]
        
        for p in self.tree_primes:
            trigger_idx = p - 2
            if 0 <= trigger_idx < max_n:
                self.check_map[trigger_idx].append(p)
        
        print(f"Optimized: Tree Primes: {len(self.tree_primes)} | Leaf Primes: {len(self.leaf_primes)}")


    def get_primes(self, limit):
        sieve = [True] * (limit + 1)
        for num in range(2, int(limit**0.5) + 1):
            if sieve[num]:
                for multiple in range(num*num, limit + 1, num):
                    sieve[multiple] = False
        return [num for num in range(2, limit + 1) if sieve[num]]


    def decode_mask(self, mask):
        components = []
        for k in range(1, self.max_n + 2):
            if (mask >> k) & 1:
                components.append(f"{k}!")
        return " + ".join(components) if components else "0"


    def search(self, start_kernel):
        print(f"\n--- FAST SEARCH STARTING ---")
        
        # Initialize Residues
        current_residues = {p: 0 for p in self.tree_primes}
        
        kernel_sum = mpz(0)
        kernel_mask = 0
        start_idx = 0
        
        # Stack: (idx, current_sum, residues_dict, weight, mask)
        stack = [(start_idx, kernel_sum, current_residues, 1.0, kernel_mask)]
        
        start_time = time.time()
        last_print = time.time()
        nodes_checked = 0
        progress = 0.0
        
        while stack:
            idx, curr_sum, residues, weight, mask = stack.pop()
            
            # --- LEAF NODE (Or End of Range) ---
            if idx >= self.max_n:
                # CHANGED: Remove + 1 for A014597 general case
                candidate = curr_sum
                
                # ADDED: Skip k=0 (empty sum)
                if candidate == 0:
                    progress += weight
                    continue
                
                possible = True
                
                # Fast Mod 64 check first
                if (candidate & 63) not in {0,1,4,9,16,17,25,33,36,41,49,57}:
                    possible = False
                
                # Check extra Leaf Primes
                if possible:
                    for p in self.leaf_primes:
                        if not self.qr_lookups[p][candidate % p]:
                            possible = False
                            break
                
                if possible:
                    if is_square(candidate):
                        root = int(candidate**0.5)
                        eqn = self.decode_mask(mask)
                        # CHANGED: Removed "+ 1" from print statement
                        print(f"\n>>> HIT: {eqn} = {candidate} ({root}^2)")
                        self.results.append((root, eqn))
                
                progress += weight
                continue


            # --- REPORTING ---
            nodes_checked += 1
            if nodes_checked & 0x3FFFF == 0:
                now = time.time()
                if now - last_print > 0.5:
                    elapsed = now - start_time
                    speed = nodes_checked / elapsed
                    eta_str = "Calc..."
                    if progress > 0:
                        eta = (elapsed / progress) - elapsed
                        eta_str = time.strftime("%H:%M:%S", time.gmtime(eta))
                    
                    print(f"\rDepth: {idx:<2} | Nodes: {nodes_checked/1000000:.1f}M | Speed: {speed/1000:.0f}k/s | Done: {progress*100:6.4f}% | ETA: {eta_str} ", end="")
                    last_print = now


            # --- PREPARE NEXT ---
            next_weight = weight * 0.5
            term_val = self.factorials[idx]
            
            primes_to_check = self.check_map[idx]
            
            # --- BRANCH 1: ADD (idx+1)! ---
            valid_add = True
            
            for p in primes_to_check:
                r = (residues[p] + self.fact_mods[idx][p]) 
                if r >= p: r -= p
                
                # CHANGED: Check residue directly (no +1)
                if not self.qr_lookups[p][r % p]:
                    valid_add = False
                    break
            
            if valid_add:
                new_residues = residues.copy()
                for p in self.tree_primes:
                    nr = new_residues[p] + self.fact_mods[idx][p]
                    if nr >= p: nr -= p
                    new_residues[p] = nr
                
                new_mask = mask | (1 << (idx+1))
                stack.append((idx + 1, curr_sum + term_val, new_residues, next_weight, new_mask))
            else:
                progress += next_weight


            # --- BRANCH 2: SKIP (idx+1)! ---
            valid_skip = True
            
            for p in primes_to_check:
                # CHANGED: Check residue directly (no +1)
                if not self.qr_lookups[p][residues[p] % p]:
                    valid_skip = False
                    break
            
            if valid_skip:
                stack.append((idx + 1, curr_sum, residues, next_weight, mask))
            else:
                progress += next_weight


        total_time = time.time() - start_time
        print(f"\n[DONE] Finished in {total_time:.2f}s. Total Nodes: {nodes_checked:,}")


    def print_summary(self):
        print("\n" + "="*60)
        print("FINAL SUMMARY (Sorted by 'a')")
        print("-" * 60)
        sorted_results = sorted(self.results, key=lambda x: x[0])
        if not sorted_results:
            print("No solutions found.")
        else:
            for root, eqn in sorted_results:
                # CHANGED: Removed "+ 1" from summary
                print(f"a = {root:<6} | {eqn} = {root}^2")
        print("="*60)


if __name__ == "__main__":
    searcher = UltraFastSearcher(MAX_N)
    for kernel in START_KERNELS:
        searcher.search(kernel)
    searcher.print_summary()
