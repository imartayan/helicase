pub trait Lexer {
    type Input;

    fn input(&self) -> &Self::Input;
}
