/* LAST UPDATED : 2022-10-16-1811 EDT */

// these are functions for numerical evaluation
let rafficot = {};

// methods with no dependencies will have filenames starting with 1000

// those with dependencies will 2000+

rafficot.create_column_vector = function(n_rows, character = 'a', indexed = true, indexed_from_zero = false) {

  let output_matrix = [];
  
  for (let y = 0; y < n_rows; y++) {
  
    output_matrix[y] = []
    
    output_matrix[y][0] = character + '<sub>' + (y+1) + '</sub>';
    if (indexed_from_zero) {
      output_matrix[y][0] = character + '<sub>' + (y) + '</sub>';
    }

    
  }
  
  return output_matrix;
};



rafficot.create_matrix = function(n_rows, n_cols) {


  let output_matrix = [];
  let arr = 'abcdefghijklmnopqrstuvwxyz';
  
  
  for (let y = 0; y < n_rows; y++) {
  
    output_matrix.push(new Array());
    
    for (let x = 0; x < n_cols; x++) {
      // output_matrix[y].push(0);
      output_matrix[y].push(arr[y*n_cols + x]);
    }
  }
  
  
  
  return output_matrix;
};



rafficot.create_matrix_ols = function(n_rows, n_cols, intercept = false) {

  let n_rows_ = n_rows;
  let n_cols_ = n_cols;
  let intercept_ = intercept;
  
  let output_matrix = [];
  let arr = 'abcdefghijklmnopqrstuvwxyz';
  
  
  
  for (let y = 0; y < n_rows; y++) {
  
    output_matrix.push(new Array());

    
    for (let x = 0; x < n_cols; x++) {
      
      let entry = 'x<sub>' + (x+1) + ',' + (y+1) + '</sub>';
      
      if (intercept_) {
        
        entry = 'x<sub>' + (x) + ',' + (y+1) + '</sub>';
        if (x === 0) {
          entry = 1;
        }

      }
      
      output_matrix[y][x] = entry;
      
    }
  }
  
  return output_matrix;
};



// to copy a matrix
rafficot.get_clone = function(matrix) {

  let arr_ = new Array();
  for (let y = 0; y < matrix.length; y++) {
    arr_[y] = [];
    for (let x = 0; x < matrix[0].length; x++) {
      arr_[y][x] = matrix[y][x]; 
    }
  }
  
  return arr_;
};


rafficot.get_col = function(matrix, x) {

 let arr = [];
 for (let y = 0; y < matrix.length; y++) {
   arr.push([matrix[y][x]]);
 }
 return arr;

};


function return_dot_product_as_number(arr_y, arr_x) {

  // the must have the same length
  let n = arr_x.length;
  if (n !== arr_y.length) {
    console.log('you cant dot product 2 vectors if they dont have the same length');
    return;
  }
  
  let str = "<p>";
  
  for (let i = 0; i < n; i++) {
    str += arr_y[i] + " &middot; " + arr_x[i];
    
    if (i !== (n-1)) {
      str += " + ";
    }
    
  }
  
  str += "</p>"
  
  return str;
};

rafficot.get_row = function(matrix, y) {

 return matrix.slice(y, y+1);

};

rafficot.get_difference = function(matrix_a__, matrix_b__) {

  let matrix_a = this.get_clone(matrix_a__);
  let matrix_b = this.get_clone(matrix_b__);

  for (let y = 0; y < matrix_a.length; y++) {
    for (let x = 0; x < matrix_a[y].length; x++) {
    
      matrix_a[y][x] -= matrix_b[y][x]; 
    }
  }
  
  return matrix_a;

};



rafficot.get_merge_left_right = function(matrix_a__, matrix_b__) {

  let matrix_a = this.get_clone(matrix_a__);
  let matrix_b = this.get_clone(matrix_b__);

  if (matrix_a.length === matrix_b.length) {
    
    for (let i = 0; i < matrix_a.length; i++) {
        
        matrix_a[i].push(matrix_b[i][0]);
        
        
    }
    
    return matrix_a;
    
  }
  
  return;

};


rafficot.get_sum = function(matrix_a__, matrix_b__) {

  let matrix_a = this.get_clone(matrix_a__);
  let matrix_b = this.get_clone(matrix_b__);

  for (let y = 0; y < matrix_a.length; y++) {
    for (let x = 0; x < matrix_a[y].length; x++) {
    
      matrix_a[y][x] += matrix_b[y][x]; 
    }
  }
  
  return matrix_a;

};


rafficot.get_sub_matrix = function(matrix, y, x) {

  // clone the matrix first
  let m = [];
  for (let i = 0; i < matrix.length; i++) {
    m.push([]);
    for (let j = 0; j < matrix[i].length; j++) {
      m[i].push(matrix[i][j]);
    }
  }

  // DELETE THE ROW
  m.splice(y, 1); // deletes the first row
  
  // DELETE THE COLUMN
  for (let i = 0; i < m.length; i++) {
    m[i].splice(x, 1);
  }
  
  return m;
};
rafficot.get_transpose = function(matrix) {

  let n_rows = matrix.length;
  let n_cols = matrix[0].length;
  
  let output_matrix = this.create_matrix(n_cols, n_rows);
  
  for (let y = 0; y < n_cols; y++) {
    for (let x = 0; x < n_rows; x++) {
    
      output_matrix[y][x] = matrix[x][y];
    
    }
  }

  return output_matrix;

};



  rafficot.get_determinant = function(matrix) {

    // exit condition 1 
    if (matrix.length === 1) {
      return matrix[0][0];
    }
  
    // exit condition 2 
    if (matrix.length === 2) {
      let a = matrix[0][0];
      let d = matrix[1][1];
      let b = matrix[0][1];
      let c = matrix[1][0];
      return a * d - b * c;
    }
    
  
    // IF THE MATRIX IS BIGGER THAN 2X2
    let y = 0;
    let det = 0;
    for (let x = 0; x < matrix[y].length; x++) {
      
      let sub_matrix = this.get_sub_matrix(matrix, y, x);  // another function I wrote, so possibly we can expand along any row-column combo
      
      let sub_det = this.get_determinant(sub_matrix);
      
      if ((sub_det + 0) === sub_det && ((matrix[y][x] + 0) === matrix[0][x])) {
        
        det += (-1)**((y+1)+(x+1)) * matrix[y][x] * this.get_determinant(sub_matrix); // is recursive 
      
      }
      
    }
    
    // show error, that if some value isnt a number, to use rafficot.show_determinant_expression()
    
    return det;

}; // closing fn
  


rafficot.get_product = function(matrix_a, matrix_b) {

  // i maybe should check that the number of columns from matrix_a matches the number of rows from matrix_b ***
  
  // the number of rows comes from matrix_a
  let n_rows = matrix_a.length;
  
  // the number of columns comes from matrix_b
  let n_cols = matrix_b[0].length;
  
  // the output matrix
  let output_matrix = this.create_matrix(n_rows, n_cols);
  
  let n = matrix_a[0].length; // ***


  // CHECK IF THE ENTIRE MATRIX IS NUMERICAL
  
  let s = 0;
  for (let y = 0; y < matrix_a.length; y++) {
    for (let x = 0; x < matrix_a[0].length; x++) { 
      s += matrix_a[y][x];   
    }   
  }
  
  for (let y = 0; y < matrix_b.length; y++) {
    for (let x = 0; x < matrix_b[0].length; x++) { 
      s += matrix_b[y][x];   
    }   
  }
  
  if ((s + 0) === s) {
    
    // loop over the rows of matrix_a
    for (let y = 0; y < n_rows; y++) {
      
      for (let x = 0; x < n_cols; x++) {
        
        output_matrix[y][x] = 0;

        // element (y,x) of the output_matrix is the dot product of row-y of matrix_a, and column-x of matrix_b
        // thats why we needed this ***
        for (let z = 0; z < n; z++) {
          output_matrix[y][x] += matrix_a[y][z] * matrix_b[z][x];
        }
   
      }
    }
    
    // console.log(matrix_a);
    // console.log(matrix_b);
    
    return output_matrix;
  }



  // IF IT IS NOT NUMERICAL, THEN...
  
    // loop over the rows of matrix_a
    for (let y = 0; y < n_rows; y++) {
      
      for (x = 0; x < n_cols; x++) {
        
        output_matrix[y][x] = '';

        // element (y,x) of the output_matrix is the dot product of row-y of matrix_a, and column-x of matrix_b
        // thats why we needed this ***
        for (let z = 0; z < n; z++) {
          
          
          if (z !== 0) {
            output_matrix[y][x] += ' + ';
          }
          
          output_matrix[y][x] += matrix_a[y][z] + '&middot;' + matrix_b[z][x];
          
          
          
          
        }
        
        
        
   
      }
    }



return output_matrix;

  

};

// RETURNS A MODIFIED COPY OF THE ORIGINAL
rafficot.get_cofactor_matrix = function(matrix__) {

  let matrix = this.get_clone(matrix__);

  if (matrix.length === 1) {
    return matrix;
  }
  
  let output_matrix = [];
  
  for (let y = 0; y < matrix.length; y++) {
    
    output_matrix.push([]);
    for (let x = 0; x < matrix[y].length; x++) {
    
      let sub_matrix = this.get_sub_matrix(matrix, y, x);
      let det = this.get_determinant(sub_matrix);
      output_matrix[y].push((-1)**((y+1)+(x+1)) * det);
      
    }
  }

  return output_matrix;

};
rafficot.get_adjugate = function(matrix__) {

  let matrix = this.get_clone(matrix__);
  
  // THE ADJUGATE OF A 1X1 MATRIX IS [[1]] OR I_{1} : THIS IS ODD TO ME
  if (matrix.length === 1) {
    return [[1]];
  }
 
  let cofactor_matrix = this.get_cofactor_matrix(matrix); 
  let adjugate_matrix = this.get_transpose(cofactor_matrix);

  return adjugate_matrix;

};
rafficot.get_inverse = function(matrix__) {


  let matrix = this.get_clone(matrix__);
  
  // IF IT IS 1X1, YOU CAN JUST RETURN THE RECIPROCAL OF THE ONE ENTRY
  if (matrix.length === 1) {
    return [[1/matrix[0][0]]];
  }
  
  // IF THERE ARE 2 OR MORE
  let adjugate_matrix = this.get_adjugate(matrix);
  let determinant = this.get_determinant(matrix);

  let n = matrix.length;  
  let matrix_inverse = this.create_matrix(n, n);
  
  for (let y = 0; y < n; y++) {
    for (let x = 0; x < n; x++) {
      matrix_inverse[y][x] = (1/determinant) * adjugate_matrix[y][x];
    }
  }
  
  return matrix_inverse;
  
};

rafficot.get_ols_coefficients = function(matrix_y__, matrix_x__) {

  let y_ = this.get_clone(matrix_y__);
  let x_ = this.get_clone(matrix_x__);
  
  let n_ = y_.length;     // number of observations
  let k_ = x_[0].length;  // number of variables including the constant
  
  let xT_ = this.get_transpose(x_);
  let xTx_ = this.get_product(xT_, x_);
  let xTx_inv_ = this.get_inverse(xTx_);
  let xTy_ = this.get_product(xT_, y_);
  
  let b_hat_ = this.get_product(xTx_inv_, xTy_);

  let yT_ = this.get_transpose(y_);
  // let yTy_ = this.get_product(yT_, y_);
  
  let y_hat_ = this.get_product(x_, b_hat_);

  let e_ = this.get_difference(y_, y_hat_); 
  let eT_ = this.get_transpose(e_);
  let eTe_ = this.get_product(eT_, e_);
  
  // let TSS_ = yTy_[0][0];
  
  // what if i combine y and x, then do the covariance matrix
  let y_bar_ = 0;
  let y_ss_ = 0;
  for (let i = 0; i < y_.length; i++) {
    y_bar_ += y_[i][0];
  }
  y_bar_ = y_bar_ / y_.length;
  
  for (let i = 0; i < y_.length; i++) {
    y_ss_ += (y_[i][0] - y_bar_)**2;
  }
  
  let TSS_ = y_ss_; // sum of squares
  let RSS_ = eTe_[0][0];

  let r2_ = 1 - RSS_ / TSS_;   // also r2 = cor(y, y_hat) = var(y_hat) / var(y)
  
  let df_ = n_ - k_;
  let residual_SE_ = Math.sqrt(RSS_ / df_);

  return {
   'b_hat':b_hat_,
   'TSS':TSS_,
   'RSS':RSS_,
   'y_bar':y_bar_,
   'n':n_,  // number of observations
   'k':k_,  // number of variables including the constant
   'df':df_,
   'r2':r2_,   // or var(y_hat) / var(y)
   'residual_SE':residual_SE_
  }
  
};




rafficot.show_matrix = function(obj) {


  // CHECK IF THE ARGUMENT IS AN ARRAY OR AN OBJECT
  
  let arr;
  
  // IF THE ARGUMENT SUPPLIED IS A 2D-ARRAY
  if (Array.isArray(arguments[0])) {

    if (rafficot.validate_matrix(arguments[0])) {
      arr = arguments[0];
    }
    
    if (!rafficot.validate_matrix(arguments[0])) {
      console.log('While you did supply an array. It wasnt a good one.');
      let div = document.createElement('div');
      div.innerHTML = 'show_matrix error : you stepped on a biscuit. The array wasnt 2d is all.'
      return div;
    }
    
  }

  // IF THE ARGUMENT IS AN OBJECT, WITH THE 2D-ARRAY AS ITS PROPERTY
  if (arguments[0].hasOwnProperty('arr')) {

    if (rafficot.validate_matrix(arguments[0].arr)) {
      arr = arguments[0].arr;
    }
    if (!rafficot.validate_matrix(arguments[0].arr)) {
      console.log('You supplied an object.');
      console.log('But not a good one.');
      let div = document.createElement('div');
      div.innerHTML = 'show_matrix error : you stepped on a biscuit.'
      return div;
    }
  }
  


  let n_rows = arr.length;
  let n_cols = arr[0].length;
  
  //console.log(obj);
  //console.log(n_rows);
  //console.log(n_cols);
  
  // THE TYPE OF BORDER : B, V, null
  
  let border_type = 'B'; // default
  if (arguments[0].hasOwnProperty('border_type')) {
    if (arguments[0].border_type === null) {
      border_type = null;
    }
    if (arguments[0].border_type === 'V') {
      border_type = 'V';
    }
  }
  
  // inner_cells is a 2d array that contains all the html td elements
  let inner_cells = [];



  let container = document.createElement('div');
  container.style.display = 'inline-block';
  container.style.verticalAlign = 'middle';

  // we need to do the border
  // and the content
  
  // EACH CELL. DEFAULT SPECS
  let border = '1px solid #f6f6f6';
  let padding = '0.5em 0.75em'; // can change
  let textAlign = "right";
  
  // START WITH THE CONTENT
  
  // MAKE THE TABLE
  let inner_table = document.createElement('table');
  inner_table.style.borderCollapse = 'collapse';
  inner_table.style.border = '0px solid transparent';
  inner_table.style.margin = 0;
  inner_table.style.padding = 0;
  inner_table.style.fontFamily = 'monospace';

  // FOR EACH ROW
  for (let y = 0; y < n_rows; y++) {
    
    inner_cells.push([]);

    let tr = document.createElement('tr');
    inner_table.appendChild(tr);

      // FOR EACH COLUMN
      for (let x = 0; x < n_cols; x++) {
        
        let td = document.createElement('td');
        td.style.padding = padding;
        td.style.border = border;
        td.style.textAlign = textAlign;
        
        inner_cells[y].push(td);
        
        tr.appendChild(td);
        
        
        
        
        // IF THE ELEMENT IS ITSELF AN ARRAY
        if (Array.isArray(arr[y][x])) {
          
          // should it inherit the other properties?
          td.appendChild(rafficot.show_matrix({'arr':arr[y][x]}));
          continue;
        }
        
        // IF THE ELEMENT IS AN HTML DOM ELEMENT
        if (arr[y][x] instanceof HTMLElement) {

          td.appendChild(arr[y][x]);
          continue;
          
        }
        
      // return obj instanceof HTMLElement;
        
        
        let cell = document.createElement('div');
        
        // IF THE ELEMENT IS AN HTML ELEMENT
        // SIMPLY APPEND IT TO CELL

          let val = arr[y][x];
          
          // IF IT IS A NUMBER
          if (val + 0 === val) {
            
            // IF WE HAVE A DECIMAL SPECIFICATION
            if (obj.hasOwnProperty('decimal_places')) {
              
              let decimal_places = obj.decimal_places;
              
              // IF THE SPECIFICATION IS AN INTEGER
              if (decimal_places%1 === 0 && decimal_places >= 0) {
                
                let k = 10**decimal_places;
                val = (Math.round(val * k) / k).toFixed(decimal_places);
                
              }
              
              
            }
            
          }
 

        cell.innerHTML = val;
        
  // HIGHLIGHTING
  if (obj.hasOwnProperty('highlight')) {
    
    if (Array.isArray(obj.highlight)) {
      
      
      let arr = obj.highlight;
      for (let i = 0; i < arr.length; i++) {
           
      if (arr[i].hasOwnProperty('col') && arr[i].hasOwnProperty('color')) {
        if (x === (arr[i].col - 1)) {
          td.style.backgroundColor = arr[i].color;
        }
      }
      
      if (arr[i].hasOwnProperty('row') && arr[i].hasOwnProperty('color')) {
        if (y === (arr[i].row - 1)) {
          td.style.backgroundColor = arr[i].color;
        }
      }      
        
      } // closing highlight-array loop
      
    } // closing if-array
    
    
    // LEGACY 
    if (obj.highlight.hasOwnProperty('row') && obj.highlight.hasOwnProperty('color')) {
      if (y === (obj.highlight.row - 1)) {
        td.style.backgroundColor = obj.highlight.color;
      }
    }
    // LEGACY ^^
    
  }     
        
        
        td.appendChild(cell);

          
        // }  
    
    } // closing x-loop
      
  } // closing y-loop
  
  // the table of values is done, now we need to add the border
  


  
  
  
  
  
  // THE INNER_TABLE (WHICH CONTAINS THE VALUES) IS NOW DONE
  
  
  
  


  /*
  
  // THIS FUNCTION RETURNS A TABLE WITH 3 ROWS, WITH 5-3-5 CELLS RESPECTIVELY
  // THE TOP AND BOTTOM ROWS HAVE 5 CELLS BECAUSE CELLS 1 AND 3 ARE THE OVERHANG FOR THE BORDERED MATRIX
  // CELL 1-1 HOLDS THE CONTENTS OF THE MATRIX 
  
  . . _ . .
  | _____ |
  . . _ . .
  
  */
  
  let border_thickness = 2; // in px
  let border_overhang = 5; // in px

  // THIS PART, THE CELLS OBJECT, IS JUST FOR EASE OF REFERENCE
  let cells = [];
  for (let y = 0; y < 3; y++) {
    cells.push([]);
    let n_cols = (y == 1) ? (3) : (5);
    for (let x = 0; x < n_cols; x++) {
      let td = document.createElement('td');
      td.style.border = 'none';
      td.style.margin = 0;
      td.style.padding = 0;
      cells[y].push(td);
    }
  }
  cells[1][1].colSpan = 3;
  
  // BORDER THICKNESS
  cells[0][0].style.height = border_thickness + 'px';
  
  
  // THE OVERHANG
  cells[0][1].style.width = '5px';
  cells[0][1].style.minWidth = '5px';
  cells[0][1].style.height = '2px';
  cells[0][1].style.minHeight = '2px';
  
  cells[0][3].style.width = '5px';
  cells[0][3].style.minWidth = '5px';
  cells[0][3].style.height = '2px';
  cells[0][3].style.minHeight = '2px';
  
  cells[2][1].style.width = '5px';
  cells[2][1].style.minWidth = '5px';
  cells[2][1].style.height = '2px';
  cells[2][1].style.minHeight = '2px';
  
  cells[2][3].style.width = '5px';
  cells[2][3].style.minWidth = '5px';
  cells[2][3].style.height = '2px';
  cells[2][3].style.minHeight = '2px';
  
  // NO BORDERS
  for (let i = 0; i < cells.length; i++) {
    for (let j = 0; j < cells[i].length; j++) {
      cells[i][j].style.backgroundColor = 'transparent';
    }
  }
    


  // THE SIDES
  cells[1][0].style.width = '2px';
  cells[1][0].style.minWidth = '2px';
  cells[1][0].style.backgroundColor = '#999';
  
  cells[1][2].style.width = '2px';
  cells[1][2].style.minWidth = '2px';
  cells[1][2].style.backgroundColor = '#999';


 
  
  if (border_type === 'B') {
    cells[0][0].style.backgroundColor = '#999'; // THE TOP LEFT OVERHANG
    cells[0][4].style.backgroundColor = '#999'; // THE TOP RIGHT OVERHANG
    
    cells[0][1].style.backgroundColor = '#999';
    cells[0][3].style.backgroundColor = '#999';
    cells[2][1].style.backgroundColor = '#999';
    cells[2][3].style.backgroundColor = '#999';
  }
  


    cells[2][0].style.backgroundColor = '#999';
    cells[2][4].style.backgroundColor = '#999';

    // NO BORDERS
    if (border_type === null) {
      for (let i = 0; i < cells.length; i++) {
        for (let j = 0; j < cells[i].length; j++) {
          cells[i][j].style.backgroundColor = 'transparent';
        }
      }
    }



  // MAKE THE ACTUAL TABLE
  let outer_table = document.createElement('table');
  outer_table.classList.add('rafficot_outer_table');
  
  outer_table.style.borderCollapse = 'collapse';  // most important, this is
  outer_table.style.border = 'none';
  outer_table.style.margin = 0;
  outer_table.style.padding = 0;

  for (let y = 0; y < cells.length; y++) {
    let tr = document.createElement('tr');
    tr.style.border = 'none';
    tr.style.margin = 0;
    tr.style.padding = 0;
    outer_table.appendChild(tr);
    
    for (let x = 0; x < cells[y].length; x++) {

      // APPLY THE CONTENTS
      if (x === 1 && y === 1) {
        cells[y][x].appendChild(inner_table);
      }
      tr.appendChild(cells[y][x]);
    }
  }


  container.appendChild(outer_table);

  // console.log(inner_cells);

  return container;

};

// THIS FILE HAS NO DEPENDENCIES, SAVE FOR R_MATRIX


// RETURNS TRUE FALSE
rafficot.validate_matrix = function(matrix) {

  
  // IF ITS NOT AN ARRAY, RETURN FALSE
  if (!Array.isArray(matrix)) {
    return false;
  };
  
  // FIRST ARRAY CHECK : CHECK THAT matrix IS AN ARRAY
  if (Array.isArray(matrix)) {
    
    // IF YES, THEN WE NEED TO CHECK THAT EACH ENTRY IS ALSO AN ARRAY
    // AND ALSO THAT THEY ARE THE SAME LENGTH
    
    let n_cols;

    for (let y = 0; y < matrix.length; y++) {

        // IF ITS NOT AN ARRAY, RETURN FALSE
        if (!Array.isArray(matrix[y])) {
          return false;
        };

        // CHECK THAT EACH ENTRY IS ALSO AN ARRAY
        if (Array.isArray(matrix[y])) {
          
          if (y === 0) {
            n_cols = matrix[0].length;  // SET n_cols
          }
          
          // CHECK THE LENGTH OF EACH ENTRY
          if (matrix[y].length !== n_cols) {
            return false;
          }
         
        } // CLOSING SECOND ARRAY CHECK

    } // CLOSING Y LOOP
    
    // IF WE CLEAR THE INNER LOOP, WE CAN RETURN TRUE
    return true;
  
  } // CLOSING FIRST ARRAY CHECK

};

rafficot.show_expression = function() {

  let arr = [[]];
  
  for (let i = 0; i < arguments.length; i++) {
    arr[0].push(arguments[i]);
  }

  return this.show_matrix({
    'arr':arr,
    'border_type':null
  });

};


rafficot.show_cofactor_matrix_expansion = function(matrix) {

  let n_rows = matrix.length;
  let n_cols = matrix[0].length;

  let cofactor_matrix = this.create_matrix(4, 4);
  
  for (let y = 0; y < n_rows; y++) {
    for (let x = 0; x < n_cols; x++) {
    
      let sign;
      if ((-1)**((y+1)+(x+1)) === 1) {
        sign = '+';
      } else {
        sign = '-'
      }
    
      let sub_matrix = this.get_sub_matrix(matrix, y, x);
     
      let minor_yx = this.show_matrix({'arr':sub_matrix,'border_type':'V'});
      // console.log(minor_yx);

      cofactor_matrix[y][x] = this.show_expression(sign, minor_yx);

    
    }
  }
  
  cofactor_matrix_el = this.show_matrix(cofactor_matrix);
  
  let arr = [['cof', this.show_matrix({arr:matrix,border_type:'B'}), '=', cofactor_matrix_el]];

  return this.show_matrix({
    'arr':arr,
    'border_type':null
  });

};

rafficot.show_determinant_expansion = function(matrix) {

  let arr = [[]];
  
  // arr[0].push(matrix);
  arr[0].push(this.show_matrix({arr:matrix,border_type:'V'}));
  arr[0].push('=');
  
  let y = 0;
  for (let x = 0; x < matrix[0].length; x++) {
    
    // let sign = (-1)**((y+1)+(x+1));
    let sign;
    if ((-1)**((y+1)+(x+1)) === 1) {
      sign = '+';
    } else {
      sign = '-'
    }
    
    if (x !== 0) {
      arr[0].push(sign);
    }
    
    arr[0].push(matrix[y][x]);
    
    let sub_matrix = this.get_sub_matrix(matrix, y, x);
    
    
    arr[0].push(this.show_matrix({arr:sub_matrix,border_type:'V'}));
    // arr[0].push('+');
    
  }

  return this.show_matrix({
    'arr':arr,
    'border_type':null
  });

};



rafficot.show_determinant_expression = function(matrix) {

    let div = document.createElement('div');
    div.style.display = 'inline-block';
    div.style.verticalAlign = 'middle';
    div.style.fontFamily = 'monospace';
      
    // exit condition 1 
    if (matrix.length === 1) {
                let div2 = document.createElement('div');
      div2.innerHTML = matrix[0][0];
      
            div.appendChild(div2);
      return div;
    }
  
    // exit condition 2 
    if (matrix.length === 2) {
      
          let div2 = document.createElement('div');
              div2.style.verticalAlign = 'middle';
      let a = matrix[0][0];
      let d = matrix[1][1];
      let b = matrix[0][1];
      let c = matrix[1][0];
      
      div2.innerHTML = a + '*' + d + ' - ' + b + '*' + c;
      div.appendChild(div2);
      
      return div;
    }
    
    
    
        let color = '#000';
        
        if ((matrix.length/3)%1 === 0) {
          color = '#58de';
        }
        
        if ((matrix.length/4)%1 === 0) {
          color = '#5c5e';
        }
        
        
    // expansion using row 1 
    let y = 0;

    for (let x = 0; x < matrix[y].length; x++) {
      
        // GET THE SIGN
        
        let div2 = document.createElement('div');
        div2.style.display = 'inline-block';
        
        let sign_el = document.createElement('div');
        let expandor_el = document.createElement('div');
        let product_operator_el = document.createElement('div');        
        
        sign_el.style.display = 'inline-block';
        
        expandor_el.style.display = 'inline-block';
        
        product_operator_el.style.display = 'inline-block';
        product_operator_el.innerHTML = '*';
        
        if (x === 0 || x !== 0) {
          
          let sign = ((-1)**((y+1)+(x+1)));
          
          if (sign === -1) {
            sign_el.innerHTML = '-';
          }
          
          if (sign === 1) {
            sign_el.innerHTML = '+'
          }

        }
        
        // NOW WE HAVE THE SIGN
        
        // console.log(matrix[y][x]);
        

        

        let parenthese = ["[", "]"]
        if (matrix.length%2 === 0) {
          parenthese = ["(", ")"]
        }
        
        let opening_bracket = document.createElement('div');
        opening_bracket.style.display = 'inline-block';
        opening_bracket.innerHTML = parenthese[0];
        opening_bracket.style.color = color;
        opening_bracket.style.fontSize = '1.' + matrix.length + 'em';
        
        let closing_bracket = document.createElement('div');
        closing_bracket.style.display = 'inline-block';
        closing_bracket.innerHTML = parenthese[1];
        closing_bracket.style.color = color;
        closing_bracket.style.fontSize = '1.' + matrix.length + 'em';
        
        let sub_matrix = this.get_sub_matrix(matrix, y, x);
        
        div.appendChild(sign_el);
        
        expandor_el.innerHTML = matrix[y][x];
        div2.appendChild(expandor_el);
        
        
        div2.appendChild(product_operator_el);
        
        div2.appendChild(opening_bracket);
        
        div2.appendChild(this.show_determinant_expression(sub_matrix));
        div2.appendChild(closing_bracket);
        
    div2.style.border = '1px solid ' + color;
    div2.style.margin = '3px';
    div2.style.padding = '3px';
      div.appendChild(div2);
        

      
    }
    

    return div;
    


    
}; // closing fn
  

rafficot.SOLVE_DET_2X2 = function(matrix) {

  let a = matrix[0][0];
  let d = matrix[1][1];
  let b = matrix[0][1];
  let c = matrix[1][0];

  let ad;
  let bc;
      
  let det;


      if ((parseFloat(a)+0) === 0 || (parseFloat(d)+0) === 0) {
        ad = 0;
      } else {
        
        // IF THEYRE BOTH NUMBERS, WE CAN MULTIPLY
        if ( (parseFloat(a)+0) === parseFloat(a) && (parseFloat(d)+0) === parseFloat(d) ) {
          ad = a * d;
        } else {
          ad = '(' + a + ')*(' + d + ')';
        }
       
      }
      
      console.log('line 91 :');
      console.log('b : ' + b);
      console.log('c : ' + c);
      if ((parseFloat(b)+0) === 0 || (parseFloat(c)+0) === 0) {
        
        bc = 0;
        
      } else {
        
        // IF THEYRE BOTH NUMBERS, WE CAN MULTIPLY
        if ( (b+0) === b && (c+0) === c ) {
          bc = b * c;
        } else {
          bc = '(' + b + ')*(' + c + ')';
        }
       
      }
      
      console.log('bc : ' + bc);
      
      if ((ad + 0) !== ad || (bc + 0) !== bc) {
        return '[' + ad + '-' + bc + ']';
      }
      
      if ((ad + 0) === ad || (bc + 0) === bc) {
        return (ad - bc); 
      }
      
} // closing fn SOLVE_DET_2X2()

rafficot.show_matrix_of_sub_matrices = function(matrix) {

  let n_rows = matrix.length;
  let n_cols = matrix[0].length;
  let output_matrix = this.create_matrix(n_rows, n_cols);

  for (let y = 0; y < matrix.length; y++) {
    
    for (let x = 0; x < matrix[y].length; x++) {
      
      let sub_matrix = this.get_sub_matrix(matrix, y, x);
      
      output_matrix[y][x] = sub_matrix;
      
    }
  }


  return rafficot.show_matrix({
    'arr':output_matrix
  });

};

rafficot.show_product = function(matrix_a, matrix_b) {

  if (!this.validate_matrix(matrix_a)) {
    return false;
  }
  
  if (!this.validate_matrix(matrix_b)) {
    return false;
  }

  // IF matrix_a HAS n COLUMNS, THEN matrix_b MUST HAVE n ROWS
  if (matrix_a[0].length !== matrix_b.length) {
    return false;
  }
  
  // THE PRODUCT
  let matrix_ab = this.get_product(matrix_a, matrix_b);

  // THE ELEMENTS
  let matrix_a_el = this.show_matrix({'arr':matrix_a});
  let matrix_b_el = this.show_matrix({'arr':matrix_b});
  let matrix_ab_el = this.show_matrix({'arr':matrix_ab});

  return this.show_expression(matrix_a_el, matrix_b_el, "=", matrix_ab_el)

};


rafficot.show_product_expansion = function(matrix_a, matrix_b) {

  if (!this.validate_matrix(matrix_a)) {
    return false;
  }
  
  if (!this.validate_matrix(matrix_b)) {
    return false;
  }
  
  // IF matrix_a HAS n COLUMNS, THEN matrix_b MUST HAVE n ROWS
  if (matrix_a[0].length !== matrix_b.length) {
    return false;
  }
  
  let n_rows = matrix_a.length;
  let n_cols = matrix_b[0].length;
  
  let matrix_ab_expansion = this.create_matrix(n_rows, n_cols);
  
  // loop over the rows of matrix_a
  for (let y = 0; y < n_rows; y++) {
    
    // loop over the columns of matrix_b
    for (let x = 0; x < n_cols; x++) {
      
      let row = this.get_row(matrix_a, y);
      let col = this.get_col(matrix_b, x);
      
      let cell = this.show_expression(this.show_matrix(row), this.show_matrix(col))
      
      matrix_ab_expansion[y][x] = cell;
      
    } // closing x-loop

  } // closing y-loop
  
  // THE ELEMENTS
  let matrix_a_el = this.show_matrix({'arr':matrix_a});
  let matrix_b_el = this.show_matrix({'arr':matrix_b});
  let matrix_ab_expansion_el = this.show_matrix({'arr':matrix_ab_expansion});

  return this.show_expression(matrix_a_el, matrix_b_el, '=', matrix_ab_expansion_el);
  
};
rafficot.show_inverse_expansion = function(matrix) {

  // let matrix_inverse = [];

  let adjugate_matrix = this.get_adjugate(matrix);
  let determinant = this.get_determinant(matrix);
  let n = matrix.length;
  
  
  let matrix_inverse = this.create_matrix(n, n);
  
  for (let y = 0; y < n; y++) {
    for (let x = 0; x < n; x++) {
      matrix_inverse[y][x] = '<div style="text-align: center;">' + adjugate_matrix[y][x] + '</div>' + '<div style="border-bottom: 1px solid black;"></div>' + '<div style="text-align: center;">' + determinant + '</div>';
    }
  }
  
  // return this.show_expression(matrix, matrix_inverse, '=', this.show_product_expansion(matrix, matrix_inverse));
  return this.show_product_expansion(matrix, matrix_inverse);

};








function DEC2EXCEL(n) {
 let num_arr = DEC2NONZEROBASE(n, 26);
 let str = '';
 num_arr.forEach(function(a) {
   str += String.fromCharCode(a+96);
 });
 return str;
};

function DEC2NONZEROBASE(n,b) {
 let arr = [];
 arr.push((n-1) % b + 1);
 if (Math.floor((n-1)/b) > 0) {
  arr = DEC2NONZEROBASE(Math.floor((n-1)/b), b).concat(arr)
 }
 return arr;
};

// egg
