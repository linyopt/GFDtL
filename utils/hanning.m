function c = hanning(m)

  if (nargin ~= 1)
    usage ('hanning (m)');
  end

  if (~ (isscalar (m) && (m == round (m)) && (m > 0)))
    error ('hanning: m has to be an integer > 0');
  end

  if (m == 1)
    c = 1;
  else
    m = m - 1;
    c = 0.5 - 0.5 * cos (2 * pi * (0 : m)' / m);
  end

end