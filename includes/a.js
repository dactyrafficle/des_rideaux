
window.addEventListener('load', function() {
  
  // TO RENDER THE KATEX
  (function() {
    
    let math = document.getElementsByClassName('math');
    for (let i = 0; i < math.length; i++) {
      katex.render(math[i].textContent, math[i]);
    }
    
  })();

  // TO COPY TEXT ON CLICK
 (function() {
   
   let arr = document.getElementsByClassName('blue');
   for (let i = 0; i < arr.length; i++) {
     arr[i].addEventListener('click', function() {
      let copyText = this.innerText;
      let input = document.createElement('input');
      input.type = 'text';
      console.log(copyText);

      input.value = copyText;
      arr[i].appendChild(input);
      input.select();
      document.execCommand('copy');
      input.remove();
     });
   }

 })();
 
});
