/*	Java reference :
		References on Java IO, Structures, etc.
*/

import java.io.*;
import java.lang.*;
import java.math.*;
import java.util.*;

/*	Regular usage:
		Slower IO :
			Scanner in = new Scanner (System.in));
			Scanner in = new Scanner (new BufferedInputStream (System.in));
			Input :
				in.nextInt () / in.nextBigInteger () / in.nextBigDecimal () / in.nextDouble ()
				in.nextLine () / in.hasNext ()
			Output :
				System.out.print (...);
				System.out.println (...);
				System.out.printf (...);
		Faster IO :
			Shown below.
		BigInteger :
			BigInteger.valueOf (int) : convert to BigInteger.
			abs / negate () / max / min / add / subtract / multiply /
				divide / remainder (BigInteger) : BigInteger algebraic.
			gcd (BigInteger) / modInverse (BigInteger mod) /
				modPow (BigInteger ex, BigInteger mod) / pow (int ex) : Number Theory.
			not () / and / or / xor (BigInteger) / shiftLeft / shiftRight (int) : Bit operation.
			compareTo (BigInteger) : comparation.
			intValue () / longValue () / toString (int radix) : converts to other types.
			isProbablePrime (int certainty) / nextProbablePrime () : checks primitive.
		BigDecimal :
			consists of a BigInteger value and a scale.
			The scale is the number of digits to the right of the decimal point.
			divide (BigDecimal) : exact divide.
			divide (BigDecimal, int scale, RoundingMode roundingMode) :
				divide with roundingMode, which may be:
					CEILING / DOWN / FLOOR / HALF_DOWN / HALF_EVEN / HALF_UP / UNNECESSARY / UP.
			BigDecimal setScale (int newScale, RoundingMode roundingMode) :
				returns a BigDecimal with newScale.
			doubleValue () / toPlainString () : converts to other types.
		Arrays :
			Arrays.sort (T [] a);
			Arrays.sort (T [] a, int fromIndex, int toIndex);
			Arrays.sort (T [] a, int fromIndex, int toIndex, Comperator <? super T> comperator);
		LinkedList <E> :
			addFirst / addLast (E) / getFirst / getLast / removeFirst / removeLast () :
				deque implementation.
			clear () / add (int, E) / remove (int) : clear, add & remove.
			size () / contains / removeFirstOccurrence / removeLastOccurrence (E) :
				deque methods.
			ListIterator <E> listIterator (int index) : returns an iterator :
				E next / previous () : accesses and iterates.
				hasNext / hasPrevious () : checks availablity.
				nextIndex / previousIndex () : returns the index of a subsequent call.
				add / set (E) / remove () : changes element.
		PriorityQueue <E> (int initcap, Comparator <? super E> comparator) :
			add (E) / clear () / iterator () / peek () / poll () / size () :
				priority queue implementations.
		TreeMap <K, V> (Comparator <? super K> comparator) :
			Map.Entry <K, V> ceilingEntry / floorEntry / higherEntry / lowerEntry (K):
				getKey / getValue () / setValue (V) : entries.
			clear () / put (K, V) / get (K) / remove (K) : basic operation.
			size () : size.
		StringBuilder :
			Mutable string.
			StringBuilder (string) : generates a builder.
			append (int, string, ...) / insert (int offset, ...) : adds objects.
			charAt (int) / setCharAt (int, char) : accesses a char.
			delete (int, int) : removes a substring.
			reverse () : reverses itself.
			length () : returns the length.
			toString () : converts to string.
		String :
			Immutable string.
			String.format (String, ...) : formats a string. i.e. sprintf.
			toLowerCase / toUpperCase () : changes the case of letters.
*/

/* Examples on Comparator :

public class Main {

	public static class Point {
		public int x;
		public int y;
		public Point () {
			x = 0;
			y = 0;
		}
		public Point (int xx, int yy) {
			x = xx;
			y = yy;
		}
	};

	public static class Cmp implements Comparator <Point> {
		public int compare (Point a, Point b) {
			if (a.x < b.x) return -1;
			if (a.x == b.x) {
				if (a.y < b.y) return -1;
				if (a.y == b.y) return 0;
			}
			return 1;
		}
	};

	public static void main (String [] args) {
		Cmp c = new Cmp ();
		TreeMap <Point, Point> t = new TreeMap <Point, Point> (c);
		return;
	}
};

*/

/*	Another way to implement is to use Comparable.
	However, equalTo and hashCode must be rewritten.
	Otherwise, containers may fail and give strange answers.
	Example :

	public static class Point implements Comparable <Point> {
		public int x;
		public int y;
		public Point () {
			x = 0;
			y = 0;
		}
		public Point (int xx, int yy) {
			x = xx;
			y = yy;
		}
		public int compareTo (Point p) {
			if (x < p.x) return -1;
			if (x == p.x) {
				if (y < p.y) return -1;
				if (y == p.y) return 0;
			}
			return 1;
		}
		public boolean equalTo (Point p) {
			return (x == p.x && y == p.y);
		}
		public int hashCode () {
			return x + y;
		}
	};
*/

//Faster IO :

public class Main {

	static class InputReader {
		public BufferedReader reader;
		public StringTokenizer tokenizer;
		public InputReader (InputStream stream) {
			reader = new BufferedReader (new InputStreamReader (stream), 32768);
			tokenizer = null;
		}
		public String next() {
			while (tokenizer == null || !tokenizer.hasMoreTokens()) {
				try {
					String line = reader.readLine();
					tokenizer = new StringTokenizer (line);
				} catch (IOException e) {
					throw new RuntimeException (e);
				}
			}
			return tokenizer.nextToken();
		}
		public BigInteger nextBigInteger() {
			return new BigInteger (next (), 10);	//	customize the radix here.
		}
		public int nextInt() {
			return Integer.parseInt (next());
		}
		public double nextDouble() {
			return Double.parseDouble (next());
		}
	}

	public static void main (String[] args) {
		InputReader in = new InputReader (System.in);

		//	Put your code here.

	}
}

