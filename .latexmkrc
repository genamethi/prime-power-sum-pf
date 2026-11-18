$pdf_mode = 1;  # Build with pdfLaTeX.
$out_dir = 'obj';
$aux_dir = $out_dir;
$pdflatex = 'pdflatex -halt-on-error -synctex=1 %O %S';
