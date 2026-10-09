:::titlepage
[[title]]
Underlying Algorithms
:::

---
# Curriculum

- $k$-mer decomposition, 
- hash functions,
- MinHash sketching. 
- Estimating the Jaccard similarity index: $$J(A, B) = \frac{\vert{}A \cap B\vert{}}{\vert{}A \cup B\vert{}} \approx \frac{\vert{}\min_s(h(A)) \cap \min_s(h(B))\vert{}}{s}$$ and converting to Mash distance: $$D \approx -\frac{1}{k} \ln\left(\frac{2J}{1+J}\right)$$

