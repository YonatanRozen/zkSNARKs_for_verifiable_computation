# zkSNARKs_for_verifiable_computation
This is a project that aims to implement a client-server model which uses verifiable computation using zkSNARKs. 

These days, computations are getting heavier due to the growing use in computationally demanding approaches
for tasks in Computer Science, such as training nerual networks for LLMs, using Hyper Performance Computing for
science needs and extraction and management of very large databases for applications. The computational power 
needed for some of these tasks can be very high and many computers cannot handle them alone, but still need the 
results for other important tasks. This means that the computation itself should be executed to some external entity
like a cloud \ remote server which *does* posses a high computational power. The computation is done remotely and the 
result is sent back to the client for future use. However, as we know, some servers are exposed to cyber weaknesses 
which means they can be hacked to, so the client cannot simply rely on the result sent back from the server. Instead, 
the server should *prove* the client the correctness of the result. In some cases the prove itself can contain sensitive
information that the server does not want to reveal, so we use a Zero Knowledge Protocol between them, which generates
a proof for the correctnes of the computation and ensures privacy of secret server data. 

A disclosure: this is my first ever experience with zkSNARKs and crypt programming in genearl, so I'll implement a 
simple system that enables the client to ask the server things like "is the equation Ax=b solveable? (and solution 
represents secret values)" or "Do you own a specific value (password, key)?". A follow up project could be to expand 
the set of possible computations to bigger problems like the training of a nerual network \ extracting information from
datasets. 
