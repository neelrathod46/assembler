"""
Genome Assembly Mutation Recovery Experiment

This script demonstrates the alignment-error correction paradox:
- Generate a random reference genome
- Introduce mutations
- Generate sequencing reads
- Apply error correction (consensus voting)
- Measure which mutations survive

Key insight: Mutations in repetitive regions (where reads appear multiple times)
are more likely to be "corrected away" by consensus-based error correction.
"""

import random
from collections import Counter, defaultdict


class GenomeAssemblyExperiment:
    def __init__(
        self,
        genome_size=10000,
        read_length=150,
        coverage_depth=10,
        mutation_rate=0.05,
        error_rate=0.005,
        # random_seed=42
    ):
        """
        Initialize the experiment parameters.
        
        Args:
            genome_size: Length of the reference genome (bp)
            read_length: Length of each read (bp)
            coverage_depth: Average coverage (number of times each base is covered)
            mutation_rate: Fraction of bases to mutate in reference (0.05 = 5%)
            error_rate: Sequencing error rate per base (0.005 = 0.5%)
            random_seed: For reproducibility
        """
        self.genome_size = genome_size
        self.read_length = read_length
        self.coverage_depth = coverage_depth
        self.mutation_rate = mutation_rate
        self.error_rate = error_rate
        self.bases = ['A', 'T', 'C', 'G']
        
        # random.seed(random_seed)
        
        # Will be populated later
        self.ref_genome = None
        self.mutated_genome = None
        self.reads_with_positions = None
        self.mutations = None
        self.consensus_genome = None
        
    def generate_reference_genome(self):
        """Generate a random reference genome."""
        print(f"\n{'='*60}")
        print("STEP 1: Generate Reference Genome")
        print(f"{'='*60}")
        
        self.ref_genome = ''.join(random.choice(self.bases) for _ in range(self.genome_size))
        print(f"Generated random genome: {self.genome_size} bp")
        print(f"First 100 bp: {self.ref_genome[:100]}")
        
        return self.ref_genome
    
    def introduce_mutations(self):
        """Introduce random mutations into the reference genome."""
        print(f"\n{'='*60}")
        print("STEP 2: Introduce Mutations")
        print(f"{'='*60}")
        
        mutated = list(self.ref_genome)
        num_mutations = int(self.genome_size * self.mutation_rate)
        
        self.mutations = []  # Track (position, original_base, new_base)
        
        for _ in range(num_mutations):
            pos = random.randint(0, self.genome_size - 1)
            original = self.ref_genome[pos]
            new_base = random.choice([b for b in self.bases if b != original])
            
            mutated[pos] = new_base
            self.mutations.append((pos, original, new_base))
        
        self.mutated_genome = ''.join(mutated)
        
        print(f"Introduced {len(self.mutations)} mutations ({self.mutation_rate*100}%)")
        print(f"Example mutations (first 5):")
        for i, (pos, orig, new) in enumerate(self.mutations[:5]):
            print(f"  Position {pos}: {orig} → {new}")
    
    def add_sequencing_errors(self, read):
        """Add random sequencing errors to a read."""
        read = list(read)
        for i in range(len(read)):
            if random.random() < self.error_rate:
                current = read[i]
                new_base = random.choice([b for b in self.bases if b != current])
                read[i] = new_base
        return ''.join(read)
    
    def generate_reads(self):
        """Generate reads from the mutated genome."""
        print(f"\n{'='*60}")
        print("STEP 3: Generate Sequencing Reads")
        print(f"{'='*60}")
        
        num_reads = (self.genome_size * self.coverage_depth) // self.read_length
        self.reads_with_positions = []
        
        for _ in range(num_reads):
            pos = random.randint(0, self.genome_size - self.read_length)
            read = self.mutated_genome[pos:pos + self.read_length]
            
            # Add sequencing errors
            read_with_errors = self.add_sequencing_errors(read)
            
            self.reads_with_positions.append((read_with_errors, pos))
        
        print(f"Generated {num_reads} reads")
        print(f"Target coverage: {self.coverage_depth}×")
        print(f"Actual average coverage: {(num_reads * self.read_length) / self.genome_size:.1f}×")
        print(f"Sequencing error rate: {self.error_rate*100}%")
        print(f"\nExample reads (first 3):")
        for i, (read, pos) in enumerate(self.reads_with_positions[:3]):
            print(f"  Read {i}: pos {pos:5d}-{pos+self.read_length:5d} | {read}")
    
    def find_duplicate_reads(self):
        """Find reads that appear multiple times (indicate repeats)."""
        print(f"\n{'='*60}")
        print("STEP 4: Find Duplicate Reads (Repetitive Regions)")
        print(f"{'='*60}")
        
        read_sequences = [r for r, p in self.reads_with_positions]
        read_counts = Counter(read_sequences)
        
        duplicates = {read: count for read, count in read_counts.items() if count > 1}
        
        print(f"Total unique reads: {len(read_counts)}")
        print(f"Reads that appear multiple times: {len(duplicates)}")
        print(f"Most common read appears {max(read_counts.values())} times")
        
        return duplicates
    
    def consensus_error_correction(self):
        """
        Apply consensus-based error correction.
        At each position, pick the most common base.
        This is what tools like error correction algorithms do.
        """
        print(f"\n{'='*60}")
        print("STEP 5: Apply Consensus Error Correction")
        print(f"{'='*60}")
        
        # Count bases at each position
        position_consensus = defaultdict(Counter)
        
        for read, pos in self.reads_with_positions:
            for i, base in enumerate(read):
                position_consensus[pos + i][base] += 1
        
        # Build consensus sequence
        self.consensus_genome = []
        for pos in range(self.genome_size):
            if position_consensus[pos]:
                # Most common base at this position
                consensus_base = position_consensus[pos].most_common(1)[0][0]
                self.consensus_genome.append(consensus_base)
            else:
                # No coverage at this position (shouldn't happen with our setup)
                self.consensus_genome.append('N')
        
        self.consensus_genome = ''.join(self.consensus_genome)
        
        print("Applied consensus voting at each position:")
        print(f"  For each base position, chose the most common base from all reads")
    
    def analyze_mutation_recovery(self):
        """
        Analyze which mutations survived error correction.
        Focus on mutations in repetitive regions.
        """
        print(f"\n{'='*60}")
        print("STEP 6: Analyze Mutation Recovery")
        print(f"{'='*60}")
        
        # For each position, count how many times it's covered
        position_coverage = defaultdict(int)
        for read, pos in self.reads_with_positions:
            for i in range(len(read)):
                position_coverage[pos + i] += 1
        
        # For each read sequence, count how many times it appears
        read_sequences = [r for r, p in self.reads_with_positions]
        read_frequency = Counter(read_sequences)
        
        # Categorize mutations
        recovered_mutations = []
        lost_mutations = []
        
        for mut_pos, original_base, mutated_base in self.mutations:
            # Check if the mutation was "corrected" or "recovered"
            if self.consensus_genome[mut_pos] == mutated_base:
                recovered_mutations.append({
                    'pos': mut_pos,
                    'original': original_base,
                    'mutated': mutated_base,
                    'coverage': position_coverage[mut_pos]
                })
            else:
                lost_mutations.append({
                    'pos': mut_pos,
                    'original': original_base,
                    'mutated': mutated_base,
                    'consensus': self.consensus_genome[mut_pos],
                    'coverage': position_coverage[mut_pos]
                })
        
        print(f"\nMutation Recovery Results:")
        print(f"  Total mutations introduced: {len(self.mutations)}")
        print(f"  Mutations recovered: {len(recovered_mutations)}")
        print(f"  Mutations lost: {len(lost_mutations)}")
        print(f"  Recovery rate: {100*len(recovered_mutations)/len(self.mutations):.1f}%")
        
        if lost_mutations:
            print(f"\nFirst 10 lost mutations (would be 'corrected away'):")
            for i, mut in enumerate(lost_mutations[:10]):
                print(f"  {i+1}. Pos {mut['pos']:5d}: "
                      f"{mut['original']}→{mut['mutated']} "
                      f"(coverage {mut['coverage']}×, "
                      f"consensus: {mut['consensus']})")
        
        return recovered_mutations, lost_mutations
    
    def analyze_repeats_vs_recovery(self):
        """
        The KEY ANALYSIS: Do mutations in repetitive regions get lost more often?
        """
        print(f"\n{'='*60}")
        print("STEP 7: Repeat Frequency vs Mutation Recovery")
        print(f"{'='*60}")
        
        # For each position, find the maximum read frequency that covers it
        position_max_repeat_freq = defaultdict(lambda: 1)
        
        read_sequences = [r for r, p in self.reads_with_positions]
        read_frequency = Counter(read_sequences)
        
        for read, pos in self.reads_with_positions:
            freq = read_frequency[read]
            for i in range(len(read)):
                position_max_repeat_freq[pos + i] = max(
                    position_max_repeat_freq[pos + i], 
                    freq
                )
        
        # Categorize mutations by repeat frequency of covering reads
        recovery_by_repeat = defaultdict(lambda: {'recovered': 0, 'lost': 0})
        
        for mut_pos, original_base, mutated_base in self.mutations:
            repeat_freq = position_max_repeat_freq[mut_pos]
            
            if self.consensus_genome[mut_pos] == mutated_base:
                recovery_by_repeat[repeat_freq]['recovered'] += 1
            else:
                recovery_by_repeat[repeat_freq]['lost'] += 1
        
        print("\nMutation recovery by repeat frequency of covering reads:")
        print(f"{'Repeat Freq':<15} {'Recovered':<12} {'Lost':<12} {'Recovery %':<12}")
        print("-" * 50)
        
        for repeat_freq in sorted(recovery_by_repeat.keys()):
            stats = recovery_by_repeat[repeat_freq]
            total = stats['recovered'] + stats['lost']
            recovery_pct = 100 * stats['recovered'] / total if total > 0 else 0
            print(f"{repeat_freq:<15} {stats['recovered']:<12} "
                  f"{stats['lost']:<12} {recovery_pct:<12.1f}%")
        
        print("\n   If recovery % DECREASES as repeat frequency increases,")
        print("   this demonstrates the alignment-error correction problem!")
        
        return recovery_by_repeat
    
    def run(self):
        """Run the complete experiment."""
        print("\n" + "="*60)
        print("GENOME ASSEMBLY MUTATION RECOVERY EXPERIMENT")
        print("="*60)
        print(f"\nParameters:")
        print(f"  Genome size: {self.genome_size} bp")
        print(f"  Read length: {self.read_length} bp")
        print(f"  Coverage depth: {self.coverage_depth}×")
        print(f"  Mutation rate: {self.mutation_rate*100}%")
        print(f"  Sequencing error rate: {self.error_rate*100}%")
        
        # Run all steps
        self.generate_reference_genome()
        self.introduce_mutations()
        self.generate_reads()
        self.find_duplicate_reads()
        self.consensus_error_correction()
        recovered, lost = self.analyze_mutation_recovery()
        repeat_analysis = self.analyze_repeats_vs_recovery()
        
        print(f"\n{'='*60}")
        print("EXPERIMENT COMPLETE")
        print(f"{'='*60}\n")
        
        return {
            'recovered': recovered,
            'lost': lost,
            'repeat_analysis': repeat_analysis
        }


if __name__ == "__main__":
    # Run the experiment with default parameters
    experiment = GenomeAssemblyExperiment(
        genome_size=10000,
        read_length=150,
        coverage_depth=10,
        mutation_rate=0.05,
        error_rate=0.005
    )
    
    results = experiment.run()
    
    # You can access results for further analysis:
    # - results['recovered']: mutations that survived
    # - results['lost']: mutations that were "corrected away"
    # - results['repeat_analysis']: breakdown by repeat frequency
