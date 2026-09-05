# Fold / Hopf Continuation

In this page, we explain how to perform continuation of Fold / Hopf points and detect the associated bifurcations.

### List of detected codim 2 bifurcation points

|Bifurcation|symbol used|
|---|---|
| Bogdanov-Takens | bt |
| Bautin | gh |
| Cusp | cusp |
| Zero-Hopf | zh |
| Hopf-Hopf | hh |

In a nutshell, all you have to do (see below) is to call `continuation(br, ind_bif, lens2)` to continue the bifurcation point stored in `br.specialpoint[ind_bif]` and set proper options.
