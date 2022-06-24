# Changelog

Here I'll record changes that I make to the signatures while consolidating `ref_element` and `HexN` into `HexRef`.
I have a feeling this will be useful when I try to get something like the `Average` example or something in Fierro to use the `HexRef` element.

* `num_dim_` is no longer set in the constructor/initialization routine; it is a static constant of 3 because the Hex element is assumed to be three-dimensional

* The members `elem_order` and `p_order` have been consolidated into `elem_order`

* `HexN::basis(...)` is now `HexRef::evaluate_basis(...)`

* `HexN::partial_xi_basis`, `HexN::partial_eta_basis`, and `HexN::partial_mu_basis` are now `HexRef::evaluate_derivative_basis_xi`, `HexRef::evaluate_derivative_basis_eta`, and `HexRef::evaluate_derivative_basis_zeta`
