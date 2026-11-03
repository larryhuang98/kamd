import { Component, type ReactNode } from "react";

interface State {
  error: Error | null;
}

export class ErrorBoundary extends Component<{ children: ReactNode; label: string }, State> {
  state: State = { error: null };

  static getDerivedStateFromError(error: Error): State {
    return { error };
  }

  render() {
    if (!this.state.error) return this.props.children;
    return (
      <div className="panel boundary" role="alert">
        <b>{this.props.label} stopped working</b>
        <span className="muted small">{this.state.error.message}</span>
        <button onClick={() => this.setState({ error: null })}>Try again</button>
      </div>
    );
  }
}
