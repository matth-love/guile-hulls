;;*-scheme-*
;;; concave-hull.scm
;;
;; Copyright (c) 2011, 2012, 2013, 2016, 2018 Matthew Love <matthew.love@colorado.edu>
;;
;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.
;;
;; The program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with the program.  If not, see <http://www.gnu.org/licenses/>.
;;
;; Commentary:
;;
;;; Code:

(define-module (hulls concave-hull)
  #:version (0 0 1)
  #:export
  (pw-concave-hull))

(define ch-pi 355/113)

(define (ch-theta point-a point-b)
  "- Scheme Procedure: ch-theta point-a point-b
    Return the angle between point-a and point-b"
  (let ((t (atan (- (vector-ref point-b 1) (vector-ref point-a 1))
		 (- (vector-ref point-b 0) (vector-ref point-a 0)))))
    (if (>= t 0) t
    	(+ t (* 2 ch-pi)))))

(define (ch-atheta point-a point-b point-z)
  "- Scheme Procedure: ch-atheta point-a point-b point-z
    Return the angle between point-z point-a and point-b"
  (let ((t (atan (- (vector-ref point-a 1) (vector-ref point-z 1))
		 (- (vector-ref point-a 0) (vector-ref point-z 0))))
	(u (atan (- (vector-ref point-b 0) (vector-ref point-a 0))
		 (- (vector-ref point-b 1) (vector-ref point-a 1)))))
    (if (>= t 0) t
    	(+ t (* 2 ch-pi)))))

(define (ch-cross o a b)
  "Three points are counter-clockwise turn if ch-cross > 0,
clockwise if ch-cross < 0 and collinear if ch-cross = 0."
  (- (* (- (car a) (car o))
	(- (cadr b) (cadr o)))
     (* (- (cadr a) (cadr o))
	(- (car b) (car o)))))

(define (ch-intersect? line-a line-b)
  "Return whether or not line-a crosses line-b.
`line-* is a pair of points: '((x y) (x y))"
  (let ((a (* (ch-cross (car line-a) (cadr line-a) (car line-b))
	      (ch-cross (car line-a) (cadr line-a) (cadr line-b))))
	(b (* (ch-cross (car line-b) (cadr line-b) (car line-a))
	      (ch-cross (car line-b) (cadr line-b) (cadr line-a)))))
    (if (and (< a 0) (< b 0)) #t #f)))

(define* (inside-polygon? point polygon #:optional (inside? #f))
  "Return #t if point is inside polygon.
'#(x y) in '(#(x1 y1 ...) #(x2 y2 ...) ... #(x1 y1 ...))
points that fall on the border are considered outside the polygon."
  (if (null? (cdr polygon)) inside?
      (let* ((p1 (car polygon)) (p2 (cadr polygon))
	     (x (vector-ref point 0)) (y (vector-ref point 1))
	     (x1 (vector-ref p1 0)) (y1 (vector-ref p1 1))
	     (x2 (vector-ref p2 0)) (y2 (vector-ref p2 1)))
	(if (not (and (>= y (min y1 y2)) (<= y (max y1 y2))
		      (<= x (max x1 x2)) (not (= y y1))))
	    (inside-polygon? point (cdr polygon) inside?)
	    (if (or (= x1 x2)
		    (<= x (+ x1 
			     (/ (* (- y y1)
				   (- x2 x1))
				(- y2 y1)))))
		(inside-polygon? point (cdr polygon) (not inside?))
		(inside-polygon? point (cdr polygon) inside?))))))

(define (pw-concave-hull points)
  "- Scheme Procedure: pw-concave-hull points distance"
  (define (<=vr1 p0 p1) 
    (<= (vector-ref p0 1) (vector-ref p1 1)))

  (if (not (pair? points)) points
      (let find-hull ((xys (sort points <=vr1))
		      (hull '()))
	(cond 
	 ((null? xys) hull)
	 ((null? hull)
	  (find-hull (cdr xys) (append (list (car xys)) hull) 0))
	 (else
	  (let ((next-point (pw-next-node (car hull) (cdr xys) last-theta (* 2 ch-pi) (car (reverse hull)))))
	    (if (null? (car next-point)) hull
		(find-hull (cdr xys) 
			   (append (list (car next-point)) hull) 
			   (cdr next-point)))))))))

(define (pw-concave-hull2 points dist)
  "- Scheme Procedure: pw-concave-hull points distance"
  (define (<=vr1 p0 p1) 
    (<= (vector-ref p0 1) (vector-ref p1 1)))
  (define (euclidean p0 p1)
    (let ((dy (- (vector-ref p1 1) (vector-ref p0 1)))
	  (dx (- (vector-ref p1 0) (vector-ref p0 0))))
      (sqrt (+ (expt dy 2)
	       (expt dx 2)))))
  (define* (pw-next-node point xys 
			 #:optional 
			 (last-theta 0) 
			 (this-theta (* 2 ch-pi)) 
			 (theta-point '()))
    (if (null? xys) (cons theta-point this-theta)
	(let* ((this-point (car xys))
	       (angle (ch-atheta point this-point theta-point))
	       (this-dist (euclidean point this-point)))
	  (cond
	   ((not (< this-dist dist)) (pw-next-node point (cdr xys) last-theta angle this-point))
	   ((and (> angle 0) (< angle this-theta))
	    (pw-next-node point (cdr xys) last-theta angle this-point))
	   (else (pw-next-node point (cdr xys) last-theta this-theta theta-point))))))
  (if (not (pair? points)) points
      (let find-hull ((xys (sort points <=vr1))
		      (hull '()) 
		      (last-theta 0))
	(cond 
	 ((null? xys) hull)
	 ((null? hull)
	  (find-hull (cdr xys) (append (list (car xys)) hull) 0))
	 (else
	  (let ((next-point (pw-next-node (car hull) (cdr xys) last-theta (* 2 ch-pi) (car (reverse hull)))))
	    (if (null? (car next-point)) hull
		(find-hull (cdr xys) 
			   (append (list (car next-point)) hull) 
			   (cdr next-point)))))))))


;; (define (closer point points close)
;;   (if (null? (cdr points)) (sort close <cadr)
;;       (let* ((this-node (car points))
;; 	     (this-dist (gds-distance (car point) (cadr point) (car this-node) (cadr this-node))))
;; 	(closer point (cdr points) (cons (list (car this-node) this-dist))))))

;; (define (expand-polygon polygon point)
;;   (let ((inside? (inside-polygon? point polygon)))
;;     ;; if point is insde the polygon, return the polygon as-is.
;;     (if inside? polygon
;; 	;; If the point is outside the polygon, add it.
;; 	;; --> find the two closest polygon nodes to the point
;; 	(let ((close-list (closer point 

;;; End
