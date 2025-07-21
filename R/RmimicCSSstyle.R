#' RmimicCSSstyle
#'
#' @description Preformatted CSS style sheet
#'
#' @author Matthew B. Pontifex, \email{pontifex@@msu.edu}, April 14, 2025
#'
#' @export
#' 
RmimicCSSstyle <- function() {
  # htmlTable::htmlTable viewer does not have ability to do accordions that I can find
  
  catout <- ''
  cssstyle <- "<style>\n" ; catout <- sprintf('%s%s', catout, cssstyle)
 
  cssstyle <- '
  :root {
    --blue: #007bff;
    --indigo: #6610f2;
    --purple: #6f42c1;
    --pink: #e83e8c;
    --red: #dc3545;
    --orange: #fd7e14;
    --yellow: #ffc107;
    --green: #28a745;
    --teal: #20c997;
    --cyan: #17a2b8;
    --white: #fff;
    --gray: #6c757d;
    --gray-dark: #343a40;
    --primary: #007bff;
    --secondary: #6c757d;
    --success: #28a745;
    --info: #17a2b8;
    --warning: #ffc107;
    --danger: #dc3545;
    --light: #f8f9fa;
    --dark: #343a40;
    --font-family-sans-serif: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, "Noto Sans", sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji";
    --font-family-monospace: SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;
  }
  
  html {
    scroll-padding-top: 75px; 
  }

  body {
    margin: 0 0 0 0;
    padding: 0 0 0 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, "Noto Sans", sans-serif, "Apple Color Emoji", "Segoe UI Emoji", "Segoe UI Symbol", "Noto Color Emoji";
    font-size: 1rem;
    font-weight: 400;
    line-height: 1.5;
    color: #000000;
    font-style: normal;
    text-align: left;
    background-color: #fff;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale; 
    width: 100%%;
    overflow-x:auto;
    overflow-y:auto;
  }
  
  hr {
    box-sizing: content-box;
    height: 0;
    overflow: visible;
  }
  
  h1, h2, h3, h4, h5, h6 {
    margin-top: 0;
    margin-bottom: 0.5rem;
  }
  
  p {
    margin-top: 0;
    margin-bottom: 1rem;
  }
  
  b,
  strong {
    font-weight: bolder;
  }
  
  sub,
  sup {
    position: relative;
    font-size: 75%%;
    line-height: 0;
    vertical-align: baseline;
  }
  
  sub {
    bottom: -.25em;
  }
  
  sup {
    top: -.5em;
  }
  
  table {
    border-collapse: collapse;
  }
  
  output {
    display: inline-block;
  }

  h1, h2, h3, h4, h5, h6,
  .h1, .h2, .h3, .h4, .h5, .h6 {
    margin-bottom: 0.5rem;
    font-weight: 500;
    line-height: 1.2;
  }
  
  h1, .h1 {
    font-size: 2.5rem;
  }
  
  h2, .h2 {
    font-size: 2rem;
  }
  
  h3, .h3 {
    font-size: 1.75rem;
  }
  
  h4, .h4 {
    font-size: 1.5rem;
  }
  
  h5, .h5 {
    font-size: 1.25rem;
  }
  
  h6, .h6 {
    font-size: 1rem;
  }
  
  hr {
    margin-top: 1rem;
    margin-bottom: 1rem;
    border: 0;
    border-top: 1px solid rgba(0, 0, 0, 0.1);
  }
  
  small,
  .small {
    font-size: 0.8em;
    font-weight: 400;
  }
  
  bigger,
  .bigger {
    font-size: 1.3em;
    font-weight: 400;
    line-height: 1.3;
  }
  
  big,
  .big {
    font-size: 1.15em;
    font-weight: 400;
    line-height: 1.15;
  }
  
  .primaryheadings {
    font-size: 1.3em;
    font-weight: 400;
    line-height: 1.3;
  } 
  .secondaryheadings {
    font-size: 0.95em;
    font-weight: 300;
    line-height: 1.15;
    color: #6c757d !important;
  } 
  }
  .basicheadinglabel {
    font-weight: 400;
    font-size: 0.6em;
    line-height: 1.15;
    color: red;
  }
  .effectlabel {
    font-weight: 400;
    font-size: 1.0em;
    line-height: 1.15;
  }
  .demographicsstext {
    font-weight: 400;
    font-size: 0.9em;
    line-height: 1.15;
  }
  .randomeffectstext {
    font-weight: 400;
    font-size: 0.9em;
    line-height: 1.15;
  }
  .fixedeffectstext {
    font-weight: 400;
    font-size: 0.9em;
    line-height: 1.15;
  }
  .explainmethodtext {
    font-weight: 300;
    font-size: 0.7em;
    line-height: 1.15;
    color: #6c757d !important;
  }
  .decompdirectiontext {
    font-weight: 500;
    font-size: 0.8em;
    line-height: 1.15;
  }
  .breakdownconstanttext {
    font-weight: 300;
    font-size: 0.9em;
    line-height: 1.15;
  }
  .posthocttesteffectstext {
    font-weight: 300;
    font-size: 0.8em;
    line-height: 1.15;
  }
  
  .breakdownapproachdiv {
    padding-top: 5px;
    padding-bottom: 0px;
  }
  .headinglabeldiv {
    padding-top: 5px;
    padding-bottom: 0px;
  }
  .effectlabeldiv {
    padding-top: 5px;
    padding-bottom: 0px;
  }
  .demographicdiv {
    padding-top: 5px;
    padding-bottom: 0px;
  }
  .anovatestresults {
    padding-top: 10px;
    padding-bottom: 10px;
  }
  .posthoctestresults {
    padding-top: 10px;
    padding-bottom: 10px;
  }
  .anovadeconstruction {
    padding-bottom: 10px;
  }
  
  .accordion-button {
     padding: 2px;
     margin: 0px;
     border: none;
     box-shadow: none;
     border-color: none;
  }         
  .accordion-button:enabled {
      background-color: white;
      color: black;
       border: none;
       box-shadow: none;
       border-color: none;
  }   
        
  .accordion-button:focus {
       border: none;
       box-shadow: none;
       border-color: transparent;
  }   
  .accordion-button:hover {
      background-color: #f8f9fa;
  }   
  .accordion-button::after {
      -webkit-filter: grayscale(1) invert(1);
      filter: grayscale(1) invert(1);
  }
  
  .text-left {
    text-align: left !important;
  }
  
  .text-right {
    text-align: right !important;
  }
  
  .text-center {
    text-align: center !important;
  }
  
  .text-justify {
    text-align: justify !important;
  }
  
  .text-muted {
    color: #6c757d !important;
  }
  
  '
  catout <- sprintf('%s%s', catout, cssstyle)
  
  cssstyle <- "</style>\n" ; catout <- sprintf('%s%s', catout, cssstyle)
  
  #cssstyle <- '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-1BmE4kWBq78iYhFldvKuhfTAU6auU8tT94WrHftjDbrCEXSU1oBoqyl2QvZ6jIW3" crossorigin="anonymous">\n' ; catout <- sprintf('%s%s', catout, cssstyle)
  
  return(catout) 
}