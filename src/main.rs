use num_bigint::BigUint;
use num_traits::Zero;
use std::io::Read;
use std::{str::FromStr, vec::Vec};
use num_complex::Complex;
use std::f64::consts::PI;
use std::error::Error;

fn represent_as_vector(value: &BigUint) -> Vec<Complex<f64>> {
  value.
    to_string().
    chars().
    rev().
    map(|symbol: char| Complex{ re: symbol.to_digit(10).unwrap() as f64, im: 0f64 }).
    collect()
}

fn calculate_length_power_of_2(length: u32) -> u32 {
  let length_of_number: u32 = (length as f64).log2().ceil() as u32;
  let power_of_2: u32 = 2u32.pow(length_of_number as u32 );

  return power_of_2
}

fn supplement(value: &mut Vec<Complex<f64>>, length: u32) {
  while (value.len() as u32) < length {
    value.push(Complex::zero());
  }
}

fn calculate_transorm_fourier(value: &Vec<Complex<f64>>, point: &Complex<f64>) -> Vec<Complex<f64>> {
  if value.len() == 1 {
    return value.clone();
  }

  let mut left = vec![Complex{ re: 0f64, im: 0f64 }; value.len() / 2];
  let mut right = vec![Complex{ re: 0f64, im: 0f64 }; value.len() / 2];

  for position in 0..(value.len() / 2) {
    left[position] = value[position * 2];
    right[position] = value[position * 2 + 1];
  }

  left = calculate_transorm_fourier(&left, &(point * point));
  right = calculate_transorm_fourier(&right, &(point * point));

  let mut power = Complex{ re: 1f64, im: 0f64 };
  let mut result= vec![Complex{ re: 0f64, im: 0f64}; value.len()];

  for position in 0..(value.len() / 2) {
    result[position as usize] = left[position] + power * right[position];
    result[position + (value.len() / 2)] = left[position] - power * right[position];
    power *= point;
  }

  result
}

fn fft_result_to_biguint(coeffs: &Vec<Complex<f64>>, base: u64) -> BigUint {
  // 1) Сначала извлекаем целые «ячейки» и делаем переносы (carry)
  let mut carry: u64 = 0;
  let mut digits: Vec<u64> = Vec::with_capacity(coeffs.len());

  for &c in coeffs {
      // округляем real-часть к ближайшему целому
      let mut v = c.re.round() as i64 + carry as i64;
      if v < 0 {
          v = 0;
      }
      let u = v as u64;
      // ячейка = остаток от деления на base
      digits.push(u % base);
      // перенос в следующий разряд
      carry = u / base;
  }
  // добираем старшие переносы
  while carry > 0 {
      digits.push(carry % base);
      carry /= base;
  }
  // удаляем ведущие нули
  while digits.len() > 1 && *digits.last().unwrap() == 0 {
      digits.pop();
  }

  // 2) Переводим в BigUint
  let mut res = BigUint::zero();
  let big_base = BigUint::from(base);
  // собираем через умножение на base и прибавление очередной ячейки
  for &dig in digits.iter().rev() {
      res = &res * &big_base + BigUint::from(dig);
  }
  res
}

fn calculate_product(left: &BigUint, right: &BigUint) -> BigUint {
  let mut left_as_vector = represent_as_vector(&left);
  let mut right_as_vector = represent_as_vector(&right);

  let length = calculate_length_power_of_2(
    (left_as_vector.len() + right_as_vector.len() - 1) as u32);
  let point = Complex::from_polar(1.0f64, 2f64 * PI / (length as f64));

  supplement(&mut left_as_vector, length);
  supplement(&mut right_as_vector, length);

  // println!("{} {}", left_as_vector.len(), right_as_vector.len());

  let vector_left_transformated = calculate_transorm_fourier(&left_as_vector, &point);
  let vector_right_transormated = calculate_transorm_fourier(&right_as_vector, &point);

  // println!("{} {}", vector_left_transformated.len(), vector_right_transormated.len());

  let mut product = vec![Complex{ re: 0f64, im: 0f64 }; vector_left_transformated.len()];

  for position in 0..vector_left_transformated.len() {
    product[position] = vector_left_transformated[position] * vector_right_transormated[position];
  }

  let point = Complex::from_polar(1.0f64, -2f64 * PI / (length as f64) );
  // вычисляем обратное преобразование
  let result_complex = calculate_transorm_fourier(&product, &point).
    into_iter().
    map(|value| {value / (length as f64)}).
    collect();

  fft_result_to_biguint(&result_complex, 10u64)
}
fn main() -> Result<(), Box<dyn Error>> {
  let mut input = String::new();

  match std::io::stdin().read_to_string(&mut input) {
    Ok(_) => {}

    Err(_) => {
      return Err("fail to read from \"stdin\"".into());
    }
  };

  let mut operands = input.split_whitespace().map(
    |input: &str | BigUint::from_str(input).expect("fail to read number")
  );

  let operand_left: BigUint = operands.next().expect("fail to get a left operand");
  let operand_right: BigUint = operands.next().expect("fail to get a right operand");
  let product: BigUint = calculate_product(&operand_left, &operand_right);

  println!("{} * {} = {}", operand_left, operand_right, product);
  Ok(())
}
