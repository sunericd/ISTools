function doOnClick(myObj) {
	if (typeof console == "object") {
		console.log("A 'click' event occured"); 
		console.log(myObj);
	}
	$(document.getElementById('success').style.display = "block")
}
/* 
function startSlides(start, end, delay) {
    setTimeout(slideshow(start,start,end, delay), delay);
}

function slideshow(frame, start, end, delay) {
    return (function() {
    $(document.getElementByClassName('.toolSlide'+ frame).fadeOut());
    if (frame == end) { frame = start; } else { frame += 1; }
    setTimeout(function(){$(document.getElementByClassName('.toolSlide'+ frame).fadeIn());}, 500);
    setTimeout(slideshow(frame, start, end, delay), delay + 500);
    })
}

startSlides(1, 4, 2000) */