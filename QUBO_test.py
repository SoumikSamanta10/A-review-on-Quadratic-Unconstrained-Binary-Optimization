from dwave.system import DWaveSampler, EmbeddingComposite
sampler_auto = EmbeddingComposite(DWaveSampler(solver={'topology__type': 'chimera'}))

linear = {('a', 'a'): -6, ('b', 'b'): -2, ('c', 'c'): -7, ('d', 'd'): -6}
quadratic = {('a', 'b'): 3, ('b', 'c'): 2, ('a', 'c'): 7, ('c', 'd'): 11}
Q = {**linear, **quadratic}

sampleset = sampler_auto.sample_qubo(Q, num_reads=50)

print(sampleset)
