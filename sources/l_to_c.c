char l_to_c(int l)
{
  switch(l) {
  case 0: return 's';
  case 1: return 'p';
  case 2: return 'd';
  case 3: return 'f';
  case 4: return 'g';
  case 5: return 'h';
  case 6: return 'i';
  default: return ( (char) l -3 + 'f');
  }
}
