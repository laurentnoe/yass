#!/usr/bin/env php
<?php

/* SVG "almost" equivalent functions (there is a need to use svg_imagecreate just followed by svg_imagefill) */
function svg_imagecreate($Xdim,$Ydim) {
    return array(
        'Xdim' => $Xdim,
        'Ydim' => $Ydim,
        'SVG' => "<svg  xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"".$Xdim."\" height=\"".$Ydim."\" viewBox=\"0 0 ".$Xdim." ".$Ydim."\" overflow=\"hidden\" "
    );
}

function svg_imagefill(&$img,$xstart,$ystart,$col) {
    $img["SVG"] .= " style=\"background-color:".$col."\">\n";
}

function svg_imagecolorallocate($img, $r_int, $g_int, $b_int) {
    return "rgb(".$r_int.",".$g_int.",".$b_int.")";
}

function svg_imageline(&$img ,$x1, $y1, $x2, $y2, $col, $width=1) {
    $x1 = round($x1, 2);
    $y1 = round($y1, 2);
    $x2 = round($x2, 2);
    $y2 = round($y2, 2);
    $img["SVG"] .= " <line x1=\"".$x1."\" y1=\"".$y1."\" x2=\"".$x2."\" y2=\"".$y2."\" style=\"stroke:".$col.";stroke-width:".$width."px;stroke-linecap:square;\"/>\n";
}

function svg_imagestring(&$img, $size, $x, $y, $str, $col) {
    $x = round($x, 2);
    $y = round($y, 2);
    $img["SVG"] .= " <text x=\"".$x."\" y=\"".$y."\" fill=\"".$col."\" text-anchor=\"left\" dominant-baseline=\"central\">".htmlspecialchars($str)."</text>\n";
}

function svg_imagestringup(&$img, $size, $x, $y, $str, $col) {
    $x = round($x, 2);
    $y = round($y, 2);
    $img["SVG"] .= " <text x=\"".$x."\" y=\"".$y."\" fill=\"".$col."\" text-anchor=\"middle\" dominant-baseline=\"central\" transform=\"rotate(90,".$x.",".$y.")\">".htmlspecialchars($str)."</text>\n";
}

function svg_imagepng(&$img,$filename) {
    $img["SVG"] .= "</svg>\n";
    file_put_contents($filename,$img["SVG"]);
}



/*
 * this is the yass (non-zoomable) dotplot for OFFLINE
 * scripting out of the yass output.
 *
 * (ugly code derived from the website php script)
 *
 *  you can of course improve it (as the rest of the software here)
 *  and propose new functionalities on github !!
 */

// GD takes more memory than such small limit
ini_set("memory_limit","-1");

// Error reporting
error_reporting(-1);
ini_set('display_errors', 'On');

if ($argc < 2) {
  die("Usage :\n".
      "  ./yass2dotplot.php  yass-output.yop\n".
      "                 { filename1=\"youpi 1\"  filename2=\"youpi 2\" }\n".
      "                 { filename=fileoutput.png }\n".
      "                 { colorForwardR=255 colorForwardG=128 colorForwardB=255 }\n".
      "                 { colorReverseR=255 colorReverseG=128 colorReverseB=255 }\n".
      "                 { colorForwardLightR=192 colorForwardLightG=224 colorForwardLightB=255 }\n".
      "                 { colorReverseLightR=255 colorReverseLightG=224 colorReverseLightB=192 }\n".
      "                 { maxAl=100000 }\n".
      "                 { svg=1 }\n".
      " => the output is by default named \"dp.png\"   (or  \"dp.svg\" if you set \"svg=1\" see before) \n".
      "    (some parameters can be changed inside this script, hacks\&fixes can also be done)\n"
      );
}
// Ugly by helpful :
// print_r($argv);

$filename = "dp.png";
$filename1="filename1";
$filename2="filename2";
$colorForwardR =  64;
$colorForwardG = 128;
$colorForwardB = 255;
$colorReverseR = 255;
$colorReverseG = 128;
$colorReverseB =  64;
$colorForwardLightR = 192;
$colorForwardLightG = 224;
$colorForwardLightB = 255;
$colorReverseLightR = 255;
$colorReverseLightG = 224;
$colorReverseLightB = 192;
$maxAl = 1000000;
$svg = 0;
/*
 * 1.1) Get the parameters
 */

for ($i = 2; $i < $argc; $i++) {
  if (!empty($argv[$i])) {
    $params = explode("\:", $argv[$i]);
    $newParams = array();
    foreach ($params as $param) {
      $nameval = explode("=", $param);
      // Ugly by helpful :
      //print_r($nameval);
      if (isset($nameval[0]) && $nameval[0] != "" && isset($nameval[1]) && $nameval[1] != "") {
        switch($nameval[0]) {
        case "filename":
          $filename = $nameval[1];
          break;
        case "filename1":
          $filename1 = $nameval[1];
          break;
        case "filename2":
          $filename2 = $nameval[1];
          break;
        case "colorForwardR":
          $colorForwardR=intval($nameval[1]);
          break;
        case "colorForwardG":
          $colorForwardG=intval($nameval[1]);
          break;
        case "colorForwardB":
          $colorForwardB=intval($nameval[1]);
          break;
        case "colorReverseR":
          $colorReverseR=intval($nameval[1]);
          break;
        case "colorReverseG":
          $colorReverseG=intval($nameval[1]);
          break;
        case "colorReverseB":
          $colorReverseB=intval($nameval[1]);
          break;
        case "colorForwardLightR":
          $colorForwardLightR=intval($nameval[1]);
          break;
        case "colorForwardLightG":
          $colorForwardLightG=intval($nameval[1]);
          break;
        case "colorForwardLightB":
          $colorForwardLightB=intval($nameval[1]);
          break;
        case "colorReverseLightR":
          $colorReverseLightR=intval($nameval[1]);
          break;
        case "colorReverseLightG":
          $colorReverseLightG=intval($nameval[1]);
          break;
        case "colorReverseLightB":
          $colorReverseLightB=intval($nameval[1]);
          break;
        case "maxAl":
          $maxAl = intval($nameval[1]);
          break;
        case "svg":
          $svg = intval($nameval[1]);
          if (($svg != 0) and ((substr($filename,-4) == '.png') or substr($filename,-4) == '.Png' or (substr($filename,-4) == '.PNG'))) {
            $filename = substr($filename, 0, -4).".svg";
          }
          break;
        }
      }
    }
  }
}


echo "[starting] \"".$argv[1]."\" -> \"".$filename."\" ... ";

/*
 * 1.0) Fix some global variables (to set the size of the dotplot)
 *      Note "dimX" can be changed to something larger
 */

$fact     = 1;
$coef     = 1000;
$espace   = 80;
$Amax     = 0;
$Bmax     = 0;
$dimX     = 750;
$dimY     = 0;
$thickness = 40;

/*
 * 2) Compute min/max coordinates to compute the display area
 */

$fp = fopen($argv[1],"r");
if (!$fp) {
  echo "";
  echo "Error : File ".($argv[1])." not Found.\n";
  return;
}

$nb = 0;
while (!feof($fp) && $nb < $maxAl) {
  $line = fgets($fp);
  /* Yass first line looks like this :
   * "*(1083406-1083463)(1290780-1290837) Ev: 9.1282 s: 58/58 f"
   */
  if (preg_match ("/\*\(([0-9]+)-([0-9]+)\)\(([0-9]+)-([0-9]+)\) Ev: .*[ ]+s:[ ]+[0-9]+\/[0-9]+[ ]+([fr])$/", $line, $matches)) {
    $a1 = intval($matches[1]);
    $a2 = intval($matches[2]);
    $b1 = intval($matches[3]);
    $b2 = intval($matches[4]);
    $format = $matches[5];

    $Amax = max($Amax,$a1);
    $Amax = max($Amax,$a2);
    $Bmax = max($Bmax,$b1);
    $Bmax = max($Bmax,$b2);
    $nb++;
  }

  /* Yass second line (if given by the -d 1, not present on -d 0) looks like this :
   * "* "first_file_name" (71473 bp) / "second file_name" (1988 bp)"
   */
  if (preg_match ("/\* \".*\" \(([0-9]+) bp\) \/ \".*\" \(([0-9]+) bp\)$/", $line, $match)) {
    $a3 = intval($match[1]);
    $b3 = intval($match[2]);
    $Amax = max($Amax, $a3);
    $Bmax = max($Bmax, $b3);
  }
}
fclose($fp);


$inverted = false;
if ($Amax < $Bmax) {
  $inverted = true;
}



/*
 * restricting the picture area to 800 pixels
 */

// Scale factor and scale coefs computations [bp/pixels]
if (!$inverted)  {
  if ($Amax == 0) {
    $Amax = 100000;
  }
  $fact = $Amax/$dimX;
  while ($coef < $Amax/10) {
    $coef *= 2;
    if ($coef < $Amax/10) {
      $coef *= 2.5;
      if ($coef < $Amax/10) {
        $coef *= 2;
      }
    }
  }
} else {
  if ($Bmax == 0) {
    $Bmax = 100000;
  }

  $fact = $Bmax/$dimX;
  while ($coef < $Bmax/10) {
    $coef *= 2;
    if ($coef < $Bmax/10) {
      $coef *= 2.5;
      if ($coef < $Bmax/10) {
        $coef *= 2;
      }
    }
  }
}

// Compute the other size (smaller) once the main one is set
if (!$inverted)  {
  $dimY = intval($Bmax/$fact);
} else {
  $dimY = intval($Amax/$fact);
}


// Picture Allocation
if ($svg == 0) {
  $image =     imagecreate(($dimX+$espace),($dimY+$espace));
} else {
  $image = svg_imagecreate(($dimX+$espace),($dimY+$espace));
}
// Colors Definition
if ($svg == 0) {
  $colBackground   =     imagecolorallocate($image, 255,255,255);
  $colForeground   =     imagecolorallocate($image,  64, 64, 64);
  $colForward      =     imagecolorallocate($image, $colorForwardR, $colorForwardG, $colorForwardB);
  $colForwardLight =     imagecolorallocate($image, $colorForwardLightR, $colorForwardLightG, $colorForwardLightB);
  $colReverse      =     imagecolorallocate($image, $colorReverseR, $colorReverseG, $colorReverseB);
  $colReverseLight =     imagecolorallocate($image, $colorReverseLightR, $colorReverseLightG, $colorReverseLightB);
} else {
  $colBackground   = svg_imagecolorallocate($image, 255,255,255);
  $colForeground   = svg_imagecolorallocate($image,  64, 64, 64);
  $colForward      = svg_imagecolorallocate($image, $colorForwardR, $colorForwardG, $colorForwardB);
  $colForwardLight = svg_imagecolorallocate($image, $colorForwardLightR, $colorForwardLightG, $colorForwardLightB);
  $colReverse      = svg_imagecolorallocate($image, $colorReverseR, $colorReverseG, $colorReverseB);
  $colReverseLight = svg_imagecolorallocate($image, $colorReverseLightR, $colorReverseLightG, $colorReverseLightB);

}
// Backgound Filing
if ($svg == 0) {
      imagefill($image, 0, 0,$colBackground);
} else {
  svg_imagefill($image, 0, 0,$colBackground);
}

// Legend Writing
if (!$inverted) {
  if ($svg == 0) {
        imagestringup($image,3,intval($dimX+60),intval($dimY/2), $filename2,$colForeground);
        imagestring  ($image,3,intval($dimX/2), intval($dimY+60),$filename1,$colForeground);
  } else {
    svg_imagestringup($image,3,intval($dimX+70),intval($dimY/2), $filename2,$colForeground);
    svg_imagestring  ($image,3,intval($dimX/2), intval($dimY+60),$filename1,$colForeground);
  }
} else {
  if ($svg == 0) {
        imagestringup($image,3,intval($dimX+60),intval($dimY/2), $filename1,$colForeground);
        imagestring  ($image,3,intval($dimX/2), intval($dimY+60),$filename2,$colForeground);
  } else {
    svg_imagestringup($image,3,intval($dimX+70),intval($dimY/2), $filename1,$colForeground);
    svg_imagestring  ($image,3,intval($dimX/2), intval($dimY+60),$filename2,$colForeground);
  }
}



/*
 * 2) Dotplot
 */


// used to memorise some "pixel based" stats to get small histograms
$x_forward = array($dimX);
$x_reverse = array($dimX);

$y_forward = array($dimY);
$y_reverse = array($dimY);

for ($x = 0; $x < $dimX; $x++) {
  $x_reverse[$x] = 0;
  $x_forward[$x] = 0;
}

for ($y = 0; $y < $dimY; $y++) {
  $y_reverse[$y] = 0;
  $y_forward[$y] = 0;
}


// dotplots alignments tracing
$nb = 0;
$fp = fopen($argv[1],"r");
while (!feof($fp) && $nb < $maxAl) {
  $line = fgets($fp);
  if (preg_match ("/\*\(([0-9]+)-([0-9]+)\)\(([0-9]+)-([0-9]+)\) Ev: .*[ ]+s:[ ]+[0-9]+\/[0-9]+[ ]+([fr])$/", $line, $matches)) {
    $a1 = intval($matches[1]);
    $a2 = intval($matches[2]);
    $b1 = intval($matches[3]);
    $b2 = intval($matches[4]);
    $format = $matches[5];

    if (!$inverted)  {
      $x1 = intval($a1/$fact);
      $x2 = intval($a2/$fact);
      $y1 = intval($b1/$fact);
      $y2 = intval($b2/$fact);
    } else {
      $y1 = intval($a1/$fact);
      $y2 = intval($a2/$fact);
      $x1 = intval($b1/$fact);
      $x2 = intval($b2/$fact);
    }

    // histograms computations
    if (preg_match("/f/",$format)) {
      if ($svg == 0) {
            imageline($image,$x1,$y1,$x2,$y2,$colForward);
      } else {
        svg_imageline($image,$x1,$y1,$x2,$y2,$colForward);
      }
      for ($x = min($x1,$x2), $y = min($y1,$y2); $x < max($x1,$x2) && $y < max($y1,$y2); $x++,$y++) {
        $x_forward[$x]++;
        $y_forward[$y]++;
      }
    } else {
      if ($svg == 0) {
            imageline($image,$x1,$y1,$x2,$y2,$colReverse);
      } else {
        svg_imageline($image,$x1,$y1,$x2,$y2,$colReverse);
      }
      for ($x = min($x1,$x2), $y = min($y1,$y2); $x < max($x1,$x2) && $y < max($y1,$y2); $x++,$y++) {
        $x_reverse[$x]++;
        $y_reverse[$y]++;
      }
    }
    $nb++;
  }
}
fclose($fp);


// max histogram computation
$xmax = 0;
for ($x = 0; $x < $dimX; $x++) {
  $xmax = max($xmax, $x_forward[$x] + $x_reverse[$x]);
}
$ymax = 0;
for ($y = 0; $y < $dimY; $y++) {
  $ymax = max($ymax, $y_forward[$y] + $y_reverse[$y]);
}

if ($xmax < $thickness && $ymax < $thickness) {
  $xmax = $thickness;
  $ymax = $thickness;
}


// Histogram plotting
for ($x = 0; $x < $dimX; $x++) {
  if (($x_forward[$x]) * $thickness/$xmax >= 1)
    if ($svg == 0) {
          imageline($image,$x,intval(1+$dimY)                                    ,$x,intval(1+$dimY+($x_forward[$x])                * $thickness/$xmax),$colForwardLight);
    } else {
      svg_imageline($image,$x,intval(1+$dimY)                                    ,$x,intval(1+$dimY+($x_forward[$x])                * $thickness/$xmax),$colForwardLight);
    }

  if (($x_reverse[$x]) * $thickness/$xmax >= 1)
    if ($svg == 0) {
          imageline($image,$x,intval(1+$dimY+($x_forward[$x]) * $thickness/$xmax),$x,intval(1+$dimY+($x_forward[$x]+$x_reverse[$x]) * $thickness/$xmax),$colReverseLight);
    } else {
      svg_imageline($image,$x,intval(1+$dimY+($x_forward[$x]) * $thickness/$xmax),$x,intval(1+$dimY+($x_forward[$x]+$x_reverse[$x]) * $thickness/$xmax),$colReverseLight);
    }
}
for ($y= 0; $y<$dimY; $y++) {
  if (($y_forward[$y]) * $thickness/$ymax >= 1)
    if ($svg == 0) {
          imageline($image,intval(1+$dimX)                                    ,$y,intval(1+$dimX+($y_forward[$y])                 * $thickness/$ymax),$y,$colForwardLight);
    } else {
      svg_imageline($image,intval(1+$dimX)                                    ,$y,intval(1+$dimX+($y_forward[$y])                 * $thickness/$ymax),$y,$colForwardLight);
    }

  if (($y_reverse[$y]) * $thickness/$ymax >= 1)
    if ($svg == 0) {
          imageline($image,intval(1+$dimX+($y_forward[$y]) * $thickness/$ymax),$y,intval(1+$dimX+($y_forward[$y]+$y_reverse[$y])  * $thickness/$ymax),$y,$colReverseLight);
    } else {
      svg_imageline($image,intval(1+$dimX+($y_forward[$y]) * $thickness/$ymax),$y,intval(1+$dimX+($y_forward[$y]+$y_reverse[$y])  * $thickness/$ymax),$y,$colReverseLight);
    }
}


// Coordinates
if (!$inverted) {
  // X coordinates
  $Xrap = $Amax/$coef;
  for ($i = 1; $i < $Xrap; $i++) {
    $reel = intval($i*$coef/$fact);
    if ($svg == 0) {
          imageline($image,$reel,$dimY,$reel,$dimY+5,$colForeground);
          imagestring($image,2,$reel,$dimY+7,$i*$coef,$colForeground);
    } else {
      svg_imageline($image,$reel,$dimY,$reel,$dimY+5,$colForeground);
      svg_imagestring($image,2,$reel,$dimY+7,$i*$coef,$colForeground);
    }
  }

  // Y coordinates
  $Yrap = $Bmax/$coef;
  for ($i = 1; $i < $Yrap; $i++) {
    $reel = intval($i*$coef/$fact);
    if ($svg == 0) {
          imageline($image,$dimX,$reel,$dimX+5,$reel,$colForeground);
          imagestring($image,2,$dimX+7,$reel,$i*$coef,$colForeground);
    } else {
      svg_imageline($image,$dimX,$reel,$dimX+5,$reel,$colForeground);
      svg_imagestring($image,2,$dimX+7,$reel,$i*$coef,$colForeground);
    }
  }
} else {
  // X coordinates
  $Xrap = $Bmax/$coef;
  for ($i = 1; $i < $Xrap; $i++) {
    $reel = intval($i*$coef/$fact);
    if ($svg == 0) {
          imageline($image,$reel,$dimY,$reel,$dimY+5,$colForeground);
          imagestring($image,2,$reel,$dimY+7,$i*$coef,$colForeground);
    } else {
      svg_imageline($image,$reel,$dimY,$reel,$dimY+5,$colForeground);
      svg_imagestring($image,2,$reel,$dimY+7,$i*$coef,$colForeground);
    }
  }

  // Y coordinates
  $Yrap = $Amax/$coef;
  for ($i = 1; $i < $Yrap; $i++) {
    $reel = intval($i*$coef/$fact);
    if ($svg == 0) {
          imageline($image,$dimX,$reel,$dimX+5,$reel,$colForeground);
          imagestring($image,2,$dimX+7,$reel,$i*$coef,$colForeground);
    } else {
      svg_imageline($image,$dimX,$reel,$dimX+5,$reel,$colForeground);
      svg_imagestring($image,2,$dimX+7,$reel,$i*$coef,$colForeground);
    }
  }
}

/* Main Frame */
if ($svg == 0) {
      imageline($image,$dimX,0,$dimX,$dimY,$colForeground);
      imageline($image,0,$dimY,$dimX,$dimY,$colForeground);
      imageline($image,0,0,0,$dimY,$colForeground);
      imageline($image,0,0,$dimX,0,$colForeground);
} else {
  svg_imageline($image,$dimX,0,$dimX,$dimY,$colForeground);
  svg_imageline($image,0,$dimY,$dimX,$dimY,$colForeground);
  svg_imageline($image,0,0,0,$dimY,$colForeground);
  svg_imageline($image,0,0,$dimX,0,$colForeground);
}

// Save Picture as a Png file
if ($svg == 0) {
      imagepng($image,$filename);
} else {
  svg_imagepng($image,$filename);
}
echo "[done]\n";
?>
