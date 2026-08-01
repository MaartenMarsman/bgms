# the panel of a decisive edge carries the wheel, not a stem

    Code
      describe_panel(panel)
    Output
      label     : intrusion-dreams
      subtitle  : evidence of presence
      evidence  : PIP > .99 | log BF = 12.4
      estimate  : median = 0.32 | 95% CI [0.24, 0.40]
      wheel     : 0.995
      wheel tags: (none)
      posterior : yes
      prior     : yes
      dots      : none
      window    : 0.00 to 0.48
      caption   : Density: the weight given inclusion. Pale share of the wheel: P(absent).

# the panel of an undecided edge splits its wheel

    Code
      describe_panel(panel)
    Output
      label     : intrusion-upset
      subtitle  : undecided
      evidence  : PIP = .60 | log BF = 0.4
      estimate  : median = 0.09 | 95% CI [0.04, 0.15]
      wheel     : 0.600
      wheel tags: (none)
      posterior : yes
      prior     : yes
      dots      : none
      window    : -0.03 to 0.21
      caption   : Density: the weight given inclusion. Pale share of the wheel: P(absent).

# a saturated edge prints the capped Bayes factor and a full wheel

    Code
      describe_panel(panel)
    Output
      label     : upset-physior
      subtitle  : evidence of presence
      evidence  : PIP > .99 | log BF > 10,000
      estimate  : median = 0.41 | 95% CI [0.35, 0.47]
      wheel     : 1.000
      wheel tags: (none)
      posterior : yes
      prior     : yes
      dots      : none
      window    : 0.00 to 0.55
      caption   : Density: the weight given inclusion. Pale share of the wheel: P(absent).

# a decisive absence with no included draw is a figure, not an error

    Code
      describe_panel(panel)
    Output
      label     : a-b
      subtitle  : evidence of absence
      evidence  : PIP < .01 | log BF = -7.2
      estimate  : (none)
      wheel     : 0.000
      wheel tags: (none)
      posterior : no
      prior     : yes
      dots      : none
      window    : -1.64 to 1.64
      caption   : No retained draw included this edge. Pale share of the wheel: P(absent).

# without edge selection the panel is the Savage-Dickey figure

    Code
      describe_panel(decisive)
    Output
      label     : intrusion-dreams
      subtitle  : no edge selection
      evidence  : log BF = 7.3
      estimate  : median = 0.32 | 95% CI [0.24, 0.40]
      wheel     : 0.999
      wheel tags: data|H1 / data|H0
      posterior : yes
      prior     : yes
      dots      : 0.399, 0.000
      window    : 0.00 to 0.49
      caption   : Savage-Dickey ratio at zero (grey dots). Accented share: P(edge | data), equal prior odds.

---

    Code
      describe_panel(absent)
    Output
      label     : intrusion-avoidth
      subtitle  : no edge selection
      evidence  : log BF = -3.5
      estimate  : median = 0.00 | 95% CI [-0.06, 0.06]
      wheel     : 0.030
      wheel tags: data|H1 / data|H0
      posterior : yes
      prior     : yes
      dots      : 0.399, 13.074
      window    : -0.14 to 0.14
      caption   : Savage-Dickey ratio at zero (grey dots). Accented share: P(edge | data), equal prior odds.

# the panel reads a Blume-Capel fit like any other

    Code
      cat("subtitle  : ", panel$subtitle, "\n", sep = "")
    Output
      subtitle  : evidence of presence
    Code
      cat("wheel tags: ", panel$wheel_labels %||% "(none)", "\n", sep = "")
    Output
      wheel tags: (none)
    Code
      cat("caption   : ", panel$caption, "\n", sep = "")
    Output
      caption   : Density: the weight given inclusion. Pale share of the wheel: P(absent).

